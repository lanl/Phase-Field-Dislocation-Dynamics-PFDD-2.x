/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

© 2022. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for
Los Alamos National Laboratory (LANL), which is operated by Triad National
Security, LLC for the U.S. Department of Energy/National Nuclear Security
Administration. All rights in the program are reserved by Triad National
Security, LLC, and the U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting on its
behalf a nonexclusive, paid-up, irrevocable worldwide license in this material 
to reproduce, prepare derivative works, distribute copies to the public, perform
 publicly and display publicly, and to permit others to do so.
------------------------------------------------------------------------- */

#include "iostream"
#include "math.h"
#include "cmath"
#include "mpi.h"
#include "string.h"
#include "complex.h"
#include "fft_fftw_slab_conduct.h"
#include "app.h"
#include "lattice.h"
#include "memory.h"
#include "error.h"
#include "dfftw_mpi.h"
#include "material.h"
#include "solve.h"

using namespace PFDD_NS;

#define DELTA(i, j)   ((i==j) ? 1 : 0)

/* ---------------------------------------------------------------------- */

FFTW_Slab_Conduct::FFTW_Slab_Conduct(PFDD_C *pfdd_p, int narg, char **arg) : FFT(pfdd_p,narg,arg)
{

  user_procgrid[0] = nprocs;
  user_procgrid[1] = user_procgrid[2] = 1;
}

/* ---------------------------------------------------------------------- */

FFTW_Slab_Conduct::~FFTW_Slab_Conduct()
{

}

/* ----------------------------------------------------------------------
initialize fftw
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::init()
{
  // Check if nx, ny, nz are set

  /*create plan and iplan for fftw*/
  stress_inc = app->stoptime;
  dimension = app->dimension;
  sigma = app->sigma;
  deltasig = app->deltasig;

  create_plan();

  /* slab decomposition */

  fftwnd_mpi_local_sizes(plan, &local_x, &local_x_start,
    &local_ny_after_trans, &local_y_start_after_trans,
    &total_local_size);

    local_y = ny;
    local_z = nz;

    // allocate arrays
    allocate();

    // Check for box extremes
    double xt=0,yt=0,zt=0;
    boxxlof = boxxhif = 0;
    boxylof = boxyhif = 0;
    boxzlof = boxzhif = 0;
    double xsize[2] = {boxxlo,boxxhi};
    double ysize[2] = {boxylo,boxyhi};
    double zsize[2] = {boxzlo,boxzhi};
    for (int i=0; i<2; i++){
      for (int j=0; j<2; j++){
        for (int k=0; k<2; k++){
          xt = (xsize[i]*four[0][0] + ysize[j]*four[0][1] + zsize[k]*four[0][2])/normfour[0];
          yt = (xsize[i]*four[1][0] + ysize[j]*four[1][1] + zsize[k]*four[1][2])/normfour[1];
          zt = (xsize[i]*four[2][0] + ysize[j]*four[2][1] + zsize[k]*four[2][2])/normfour[2];
          if(xt > boxxhif) boxxhif = xt;
          if(xt < boxxlof) boxxlof = xt;
          if(yt > boxyhif) boxyhif = yt;
          if(yt < boxylof) boxylof = yt;
          if(zt > boxzhif) boxzhif = zt;
          if(zt < boxzlof) boxzlof = zt;
        }
      }
    }
    cntf[0] = (boxxhif + boxxlof)/2.0;
    cntf[1] = (boxyhif + boxylof)/2.0;
    cntf[2] = (boxzhif + boxzlof)/2.0;
    xprdf = fabs(boxxhif - boxxlof);
    yprdf = fabs(boxyhif - boxylof);
    zprdf = fabs(boxzhif - boxzlof);

  material->dslip = app->sc * material->dslip;
  if (fft->me == 0) {
    if (screen) fprintf(screen,"dslip = %lf\n", material->dslip);
    if (logfile) fprintf(logfile,"dslip = %lf\n", material->dslip);
  }
  }

  /* ----------------------------------------------------------------------
  setup fftw
  ------------------------------------------------------------------------- */

  void FFTW_Slab_Conduct::setup()
  {
    rotate_stress();
    frec();
    stiffness();
    if(app->tflag == 1) stiffness_CTB();
    greens_function();
    Bmatrix();
    // Fmatrix();
  }
  /* ----------------------------------------------------------------------
  setup global box
  assumes boxlo/hi are already set
  ------------------------------------------------------------------------- */

  void FFTW_Slab_Conduct::create_plan()
  {

    /*create plan and iplan for fftw*/

    plan = fftw3d_mpi_create_plan(world, nx, ny, nz, FFTW_FORWARD, FFTW_ESTIMATE);

    iplan = fftw3d_mpi_create_plan(world, nx, ny, nz, FFTW_BACKWARD, FFTW_ESTIMATE);

  }

  /* ----------------------------------------------------------------------
  allocate
  ------------------------------------------------------------------------- */

  void FFTW_Slab_Conduct::allocate()
  {

    slip_systems = app->slip_systems;
    num_planes = app->num_planes;
    dimension = app->dimension;
    int npa = app->nphase-1;
    memory->create(conddata_fftw,total_local_size,"conddata_fftw");
    memory->create(conddata_real,total_local_size,"conddata_real");
    memory->create(conddata_bpart,total_local_size,"conddata_bpart");
    memory->create(condtemp_data,total_local_size,"condtemp_data");
    memory->create(jcurrent,total_local_size,"jcurrent");
    memory->create(jcurrent_fftw,total_local_size,"jcurrent_fftw");
    memory->create(javerage,dimension,"javerage");
    memory->create(data_fftw,total_local_size*slip_systems,"data_fftw");
    memory->create(data_real,total_local_size*slip_systems,"data_real");
    memory->create(data_real_2,total_local_size*slip_systems,"data_real_2");
    memory->create(temp_data,total_local_size*slip_systems,"temp_data");
    memory->create(work,total_local_size,"work");
    memory->create(work_strain,total_local_size,"work_strain");
    memory->create(data_core,total_local_size*slip_systems,"data_core");
    memory->create(data_strain,total_local_size*dimension*dimension,"data_strain");
    memory->create(data_strain_NoExt,total_local_size*dimension*dimension,"data_strain_NoExt");

    memory->create(AA,total_local_size*dimension*dimension*dimension*dimension,"AA");
    if (app->nphase > 1) {
      memory->create(vstrain_complex,total_local_size*dimension*dimension,"vstrain_complex");
      memory->create(elastic_energy_den_vs,total_local_size*dimension*dimension,"elastic_energy_den_vs");
      memory->create(phase,local_x*local_y*local_z,"phase");
      memory->create(vstrain,2.0*dimension*dimension*local_x*local_y*local_z,"vstrain");
      memory->create(MM,npa,dimension,dimension,dimension,dimension,"MM");
      memory->create(EE,total_local_size*dimension*dimension*dimension*dimension,"EE");
      memory->create(ave_Vs,dimension,dimension,"ave_Vs");
      memory->create(aveVs,dimension,dimension,"aveVs");
      for (int i=0; i<local_x*local_y*local_z*dimension*dimension; i++){
        vstrain_complex[i].re = 0.0;
        vstrain_complex[i].im = 0.0;
      }
      for (int i=0; i<2*local_x*local_y*local_z*dimension*dimension; i++){
        vstrain[i] = 0.0;
      }
    }

    memory->create(xi,norder,2.0*slip_systems*local_x*local_y*local_z,"xi");
    memory->create(condxi,condnorder,2*local_x*local_y*local_z,"condxi"); // 1 scalar order parameter for charge density
    memory->create(xi_sum,norder,2.0*slip_systems*local_x*local_y*local_z,"xi_sum");
    memory->create(xi_temp,norder,2.0*slip_systems*local_x*local_y*local_z,"xi_temp");
    for(int i=0; i<norder; i++){
      for(int j=0; j<2.0*slip_systems*local_x*local_y*local_z; j++){
      xi[i][j] = 0;
      xi_temp[i][j] = 0;
      }
      for(int j=0; j<2.0*local_x*local_y*local_z; j++)
      xi_sum[i][j] = 0.0;
    }
    memory->create(dreg,norder,2.0*num_planes*local_x*local_y*local_z,"dreg");
    memory->create(dreg_sum,norder,2.0*local_x*local_y*local_z,"dreg_sum");
    for(int i=0; i<norder; i++){
      for(int j=0; j<2.0*num_planes*local_x*local_y*local_z; j++)
        dreg[i][j] = 0;
      for(int j=0; j<2.0*local_x*local_y*local_z; j++)
        dreg_sum[i][j] = 0.0;
    }
    memory->create(dregp,norder,2.0*num_planes*local_x*local_y*local_z,"dregp");
    memory->create(dregp_sum,norder,2.0*local_x*local_y*local_z,"dregp_sum");
    for(int i=0; i<norder; i++){
      for(int j=0; j<2.0*num_planes*local_x*local_y*local_z; j++)
        dregp[i][j] = 0;
      for(int j=0; j<2.0*local_x*local_y*local_z; j++)
        dregp_sum[i][j] = 0.0;
    for(int i=0; i<condnorder; i++){
      for(int j=0; j<local_x*local_y*local_z; j++)
        condxi[i][j] = 0;
      }
    }
    memory->create(xo,slip_systems*local_x*local_y*local_z,"xo");
    memory->create(xcount,local_x*local_y*local_z,"xcount");
    memory->create(fx,local_x*local_y*local_z,"fx");
    memory->create(fy,local_x*local_y*local_z,"fy");
    memory->create(fz,local_x*local_y*local_z,"fz");
    memory->create(f,slip_systems*local_x*local_y*local_z,"f");
    memory->create(r,dimension,"r");
    memory->create(C,dimension,dimension,dimension,dimension,"C");
    memory->create(S,dimension,dimension,dimension,dimension,"S");
    memory->create(C_CTB,dimension,dimension,dimension,dimension,"C_CTB");
    memory->create(S_CTB,dimension,dimension,dimension,dimension,"S_CTB");
    //memory->create(G,local_x,local_y,local_z,dimension,dimension,"G");

    memory->create(BB,slip_systems*slip_systems*local_x*local_y*local_z,"BB");
    memory->create(FF,slip_systems*total_local_size*dimension*dimension,"FF");
    memory->create(DD,slip_systems*total_local_size*dimension*dimension,"DD");

    //Grad allocation
    memory->create(gradx,slip_systems*local_x*local_y*local_z,"gradx");
    memory->create(grady,slip_systems*local_x*local_y*local_z,"grady");
    memory->create(gradz,slip_systems*local_x*local_y*local_z,"gradz");
    //theta -- angle between dislocation line tangent and Burges vector
    memory->create(theta,slip_systems*local_x*local_y*local_z,"theta");
    memory->create(hess,dimension*(dimension-1)*slip_systems*local_x*local_y*local_z,"hess");
    memory->create(nhn,slip_systems*local_x*local_y*local_z,"nhn");
    memory->create(picked,slip_systems*local_x*local_y*local_z,"picked");
    //Grad allocation FCC
    memory->create(gradx_d,num_planes*local_x*local_y*local_z,"gradx_d");
    memory->create(grady_d,num_planes*local_x*local_y*local_z,"grady_d");
    memory->create(gradz_d,num_planes*local_x*local_y*local_z,"gradz_d");
    //theta -- angle between dislocation line tangent and Burges vector
    memory->create(theta_d,num_planes*local_x*local_y*local_z,"theta_d");
    memory->create(hess_d,dimension*(dimension-1)*num_planes*local_x*local_y*local_z,"hess_d");
    memory->create(nhn_d,num_planes*local_x*local_y*local_z,"nhn_d");
    memory->create(picked_d,num_planes*local_x*local_y*local_z,"picked_d");
    memory->create(dt_kmc,num_planes*local_x*local_y*local_z,"dt_kmc");
    memory->create(firstpass,local_y,2,"firstpass");

    memory->create(fcore,num_planes*local_x*local_y*local_z,"fcore");//Check what NP is!
    memory->create(f_core,num_planes*local_x*local_y*local_z,"f_core"); // added for MPI in extended core energy by Claire 7/31/18
    memory->create(df1core,num_planes*local_x*local_y*local_z,"df1core");
    memory->create(df2core,num_planes*local_x*local_y*local_z,"df2core");
    memory->create(df3core,num_planes*local_x*local_y*local_z,"df3core");
    memory->create(dE_core,slip_systems*local_x*local_y*local_z,"dE_core");

    memory->create(E_int,local_x*local_y*local_z,"E_int");
    memory->create(E_ext,local_x*local_y*local_z,"E_ext");
    memory->create(E_lat,local_x*local_y*local_z,"E_lat");

    memory->create(xn,slip_systems,dimension,"xn");
    memory->create(xb,slip_systems,dimension,"xb");

    memory->create(data_sigma,2.0*dimension*dimension*local_x*local_y*local_z,"data_sigma");
    memory->create(data_sigma_NoExt,2.0*dimension*dimension*local_x*local_y*local_z,"data_sigma_NoExt");
    memory->create(data_asigma,2.0*dimension*dimension*local_x*local_y*local_z,"data_asigma");
    memory->create(data_tsigma,2.0*dimension*dimension*local_x*local_y*local_z,"data_tsigma");
    memory->create(data_psigma,2.0*dimension*dimension*local_x*local_y*local_z,"data_psigma");
    memory->create(data_usigma,2.0*dimension*dimension*local_x*local_y*local_z,"data_usigma");
    memory->create(data_rsigma,2.0*dimension*dimension*local_x*local_y*local_z,"data_rsigma");
    memory->create(tau,slip_systems,"tau");
    memory->create(tau_2,slip_systems*local_x*local_y*local_z,"tau_2");
    //memory->create(sigma,dimension,dimension,"sigma");
    //memory->create(deltasig,dimension,dimension,"deltasig");
    memory->create(sigma_rot,dimension,dimension,"sigma_rot");
    memory->create(local_sigma,dimension,dimension,"local_sigma");
    memory->create(ave_sigma,dimension,dimension,"ave_sigma");
    memory->create(local_sigma_NoExt,dimension,dimension,"local_sigma_NoExt");
    memory->create(ave_sigma_NoExt,dimension,dimension,"ave_sigma_NoExt");

    memory->create(eps,slip_systems,dimension,dimension,"eps");
    memory->create(data_eps,2*dimension*dimension*local_x*local_y*local_z,"data_eps");
    memory->create(data_eps_NoExt,2*dimension*dimension*local_x*local_y*local_z,"data_eps_NoExt");
    memory->create(data_epsd,2*dimension*dimension*local_x*local_y*local_z,"data_epsd");
    memory->create(avepsd,dimension,dimension,"avepsd");
    memory->create(ave_epsd,dimension,dimension,"ave_epsd");
    memory->create(ave_eps,dimension,dimension,"ave_eps");
    memory->create(avepst,dimension,dimension,"avepst");
    memory->create(avepsts,dimension,dimension,"avepsts");
    memory->create(aveps,dimension,dimension,"aveps");
    memory->create(ave_eps_NoExt,dimension,dimension,"ave_eps_NoExt");
    memory->create(aveps_NoExt,dimension,dimension,"aveps_NoExt");

    memory->create(delta,num_planes*local_x*local_y*local_z,"delta");
    memory->create(ddelta,num_planes*local_x*local_y*local_z,"ddelta");

    intave = 0;
    latave = 0;
    extave = 0;
    for(int i=0; i<dimension; i++){
      for(int j=0; j<dimension; j++){
        local_sigma[i][j] = 0;
        ave_sigma[i][j] = 0;
        local_sigma_NoExt[i][j] = 0;
        ave_sigma_NoExt[i][j] = 0;
        avepsd[i][j] = 0;
        ave_epsd[i][j] = 0;
        ave_eps[i][j] = 0;
        avepst[i][j] = 0;
        avepsts[i][j] = 0;
        aveps[i][j] = 0;
        ave_eps_NoExt[i][j] = 0;
        aveps_NoExt[i][j] = 0;
        if (app->nphase > 1) {
          ave_Vs[i][j] = 0;
          aveVs[i][j] = 0;
        }
      }
    }
    int lN1 = local_x;
    int N2 = local_y;
    int N3 = local_z;
    int NS = slip_systems;
    int NP = num_planes;
    int nmpi = NS*N2*N3;
    for (int i=0; i<nmpi*lN1; i++) {
      picked[i] = 0.0;
    }
    int nmpi_d = NP*N2*N3;
    for (int i=0; i<nmpi_d*lN1; i++) {
      picked_d[i] = 0.0;
      dt_kmc[i] = 0.0;
    }
      for(int j=0;j<N2;j++){
      for(int k=0;k<2;k++){
        firstpass[j][k]=0;
      }}
    for(int nmpi_pt=0;nmpi_pt<lN1*N2*N3;nmpi_pt++){
      xcount[nmpi_pt]=0;
    }
  }
  /* ----------------------------------------------------------------------
  Calculate sigma_rot
  ------------------------------------------------------------------------- */

  void FFTW_Slab_Conduct::rotate_stress()
  {
    if (fft->me == 0) {
      if (screen) fprintf(screen,"Rotate stress ...\n");
      if (logfile) fprintf(logfile,"Rotate stress ...\n");
    }
    int is, i, j, k, l;
    int ND = dimension;
    int ortho = 0; // 0 -- non-ortho, 1 -- ortho
    // norvol = det(prim), abs(norvol) = 1 if ortho
    if((fabs(norvol)<1.01) && (fabs(norvol)>0.99)) ortho =1;
    for(i=0; i<ND; i++){
      for(j=0; j<ND; j++){
        sigma_rot[i][j] = 0.0;
        if (ortho){
          //orthogonal grid, input stress is in local coord, need rotation
          // for(k=0; k<ND; k++)
          // for(l=0; l<ND; l++){
          //   sigma_rot[i][j] += prim[i][k]*sigma[k][l]*invprim[l][j];
          // }
          sigma_rot[i][j] = sigma[i][j];
        }
        else{
          //non-orthogonal grid, input stress is in global coord, no rotation needed
          sigma_rot[i][j] = sigma[i][j];
        }
        if (fft->me == 0) {
          if (screen) fprintf(screen,"sigam[%d][%d] = %lf,  sigma_rot[%d][%d] = %lf\n", i, j, sigma[i][j], i, j, sigma_rot[i][j]);
          if (logfile) fprintf(logfile,"sigam[%d][%d] = %lf,  sigma_rot[%d][%d] = %lf\n", i, j, sigma[i][j], i, j, sigma_rot[i][j]);
        }
      }
    }

  }

  /* ----------------------------------------------------------------------
  assign frequencies
  ------------------------------------------------------------------------- */

  void FFTW_Slab_Conduct::frec()
  {
    int i,j,k,ksym, nf;
    double d1,d2,d3;
    double kx=0.0, ky=0.0, kz=0.0;
    d1 = d2 = d3 = 1.0;
    double kxnorm = 0.0;

    for(i=0; i<local_x; i++){
      for(j=0; j<ny; j++){
        for(k=0; k<nz; k++){
          nf = k+(j)*nz+(i)*nz*ny;
          kx = local_x_start+i;
          ky = j;
          kz = k;
          /* frecuency in x */
          // if (kx >= nx/2) kx = kx - nx;
          if (kx >= nx/2) {
            kx = static_cast<double>(kx - nx)/static_cast<double>(nx);
          }
          else {
            kx = static_cast<double>(kx)/static_cast<double>(nx);
          }


          // frequency in y
          // if (ky >= ny/2) ky = ky - ny;
          if (ky >= ny/2) {
            ky = static_cast<double>(ky - ny)/static_cast<double>(ny);
          }
          else {
            ky = static_cast<double>(ky)/static_cast<double>(ny);
          }


          // frequency in z
          // if (kz >= nz/2) kz = kz - nz;
          if (kz >= nz/2) {
            kz = static_cast<double>(kz - nz)/static_cast<double>(nz);
          }
          else {
            kz = static_cast<double>(kz)/static_cast<double>(nz);
          }

          fx[nf] = (kx*four[0][0] + ky*four[1][0] + kz*four[2][0]);
          fy[nf] = (kx*four[0][1] + ky*four[1][1] + kz*four[2][1]);
          fz[nf] = (kx*four[0][2] + ky*four[1][2] + kz*four[2][2]);
          /*
          kxnorm = sqrt(fx[nf]*fx[nf] + fy[nf]*fy[nf] + fz[nf]*fz[nf]);
          if(kxnorm > 0){
          fx[nf] /= kxnorm;
          fy[nf] /= kxnorm;
          fz[nf] /= kxnorm;
        }
        */
      }
    }
  }
  /*
  int n = 0, ix = 0;
  for(j=0; j<local_x; j++){
  n = j*nz*ny;
  ix = local_x_start+j;
  if(me == 0){
  if (logfile) fprintf(logfile,"fx[%d] = %lf\n",ix, fx[n]);
  if (screen) fprintf(screen,"fx[%d] = %lf\n",ix, fx[n]);
}
if(me == 1){
if (logfile) fprintf(logfile,"fx[%d] = %lf\n",ix, fx[n]);
if (screen) fprintf(screen,"fx[%d] = %lf\n",ix, fx[n]);
}
}
for(j=0; j<ny; j++){
n = j*nz;
if(me == 0){
if (logfile) fprintf(logfile,"fy[%d] = %lf\n",j, fy[n]);
if (screen) fprintf(screen,"fy[%d] = %lf\n",j, fy[n]);
}
}
for(j=0; j<nz; j++){
n = j;
if(me == 0){
if (logfile) fprintf(logfile,"fz[%d] = %lf\n",j, fz[n]);
if (screen) fprintf(screen,"fz[%d] = %lf\n",j, fz[n]);
}
}
*/
}

/* ----------------------------------------------------------------------
Calculate stiffness
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::stiffness()
{
  double C11, C12, C44, mu,ll,young,xnu;
  double CI[6][6], SI[6][6], bb[6][6], fac[6][6], ed1[3], ed2[3], ed3[3], det_temp;
  int i, j, k, m, p, q, n,l;
  C44 = material->C44;
  C12 = material->C12;
  C11 = material->C11;
  mu  = material->mu;
  ll = material->ll;
  // double S44 = material->S44;
  // double S12 = material->S12;
  young  = material->young;
  xnu = material->xnu;
  int aflag = material->aflag;

  /* set Cijkl*/
  if (aflag == 0) {
    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        CI[i][j] = 0.;
        KCK[i][j] = 0.;
      }
    }
    CI[0][0] = C11;
    CI[0][1] = C12;
    CI[3][3] = C44;
    CI[1][1] = CI[0][0];
    CI[2][2] = CI[0][0];
    CI[0][2] = CI[0][1];
    CI[1][2] = CI[0][1];
    CI[4][4] = CI[3][3];
    CI[5][5] = CI[3][3];
    for (i=0; i<6; i++) {
       for (j=0; j<6; j++) {
         if(fabs(CI[i][j]) > 1.E-5){ //JPC
           CI[j][i] = CI[i][j];
         }
       }
     }

     // no rotation
     for (i=0; i<6; i++) {
       for (j=0; j<6; j++) {
           KCK[i][j] = CI[i][j];
       }
     }

    // BCC grid
    //  ed3[0] = -1.0/sqrt(2.0); ed3[1] =  1.0/sqrt(2.0); ed3[2] = 0.0;
    //  ed2[0] = 1.0/sqrt(6.0); ed2[1] = 1.0/sqrt(6.0); ed2[2] = -2.0/sqrt(6.0);
    //  ed1[0] = 1.0/sqrt(3.0); ed1[1] =  1.0/sqrt(3.0); ed1[2] = 1.0/sqrt(3.0);
    ed1[0] = prim[0][0]; ed1[1] =  prim[1][0]; ed1[2] = prim[2][0];
    ed2[0] = prim[0][1]; ed2[1] =  prim[1][1]; ed2[2] = prim[2][1];
    ed3[0] = prim[0][2]; ed3[1] =  prim[1][2]; ed3[2] = prim[2][2];
    rotate_stiffness(CI, ed1, ed2, ed3); 


    // GB stiffness
    // KCK[0][0]=214.30; KCK[0][1]=98.35; KCK[0][2]=89.87; KCK[0][3]=1.97; KCK[0][4]=-0.53; KCK[0][5]=-3.45;
    // KCK[1][1]=214.82; KCK[1][2]=126.25; KCK[1][3]=-0.38; KCK[1][4]=-0.13; KCK[1][5]=-23.45;
    // KCK[2][2]=210.26; KCK[2][3]=0.39; KCK[2][4]=0.52; KCK[2][5]=-11.69;
    // KCK[3][3]=67.97; KCK[3][4]=-0.04; KCK[3][5]=-0.75;
    // KCK[4][4]=31.60; KCK[4][5]=0.0;
    // KCK[5][5]=44.81;
    // for (i=0; i<6; i++) {
    //   for (j=0; j<6; j++) {
    //       KCK[j][i] = KCK[i][j];
    //   }
    // }
    // for (i=0; i<6; i++) {
    //   for (j=0; j<6; j++) {
    //       KCK[i][j] = KCK[i][j]*1.E9;
    //   }
    // }
    // GB stiffness
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {
            C[i][j][k][l] = 0.;
            S[i][j][k][l] = 0.;
          }
        }
      }
    }

    for (i = 0; i < 3; i++){
      C[0][0][i][i] = KCK[0][i];
      C[1][1][i][i] = KCK[1][i];
      C[2][2][i][i] = KCK[2][i];
      C[1][2][i][i] = KCK[3][i];
      C[2][0][i][i] = KCK[4][i];
      C[0][1][i][i] = KCK[5][i];
    }

    C[0][0][1][2] = KCK[0][3];
    C[1][1][1][2] = KCK[1][3];
    C[2][2][1][2] = KCK[2][3];
    C[1][2][1][2] = KCK[3][3];
    C[2][0][1][2] = KCK[4][3];
    C[0][1][1][2] = KCK[5][3];

    C[0][0][2][0] = KCK[0][4];
    C[1][1][2][0] = KCK[1][4];
    C[2][2][2][0] = KCK[2][4];
    C[1][2][2][0] = KCK[3][4];
    C[2][0][2][0] = KCK[4][4];
    C[0][1][2][0] = KCK[5][4];

    C[0][0][0][1] = KCK[0][5];
    C[1][1][0][1] = KCK[1][5];
    C[2][2][0][1] = KCK[2][5];
    C[1][2][0][1] = KCK[3][5];
    C[2][0][0][1] = KCK[4][5];
    C[0][1][0][1] = KCK[5][5];

    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {
            if(fabs(C[i][j][k][l]) > 1.E-5){
              C[i][j][l][k] = C[i][j][k][l];
              C[j][i][k][l] = C[i][j][k][l];
              C[j][i][l][k] = C[i][j][k][l];
            }
          }
        }
      }
    }
    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        SI[i][j] = 0.;
      }
    }
    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        bb[i][j] = 0.;
        fac[i][j] = 0.;
      }
    }
    k = 6;

    det_temp = determinant(KCK, k);
    if(me == 0){
      printf("KCK det= %E\n", det_temp);
    }
    for (q = 0;q < k; q++) {
      for (p = 0;p < k; p++) {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++) {
          for (j = 0;j < k; j++){
            if (i != q && j != p){
              bb[m][n] = KCK[i][j];
              if (n < (k - 2)) {
                n++;
              }
              else {
                n = 0;
                m++;
              }
            }
          }
        }
        fac[q][p] = pow(-1, q + p) * determinant(bb, k - 1);
      }
    }

    for (i = 0; i < k; i++){
      for (j = 0; j < k; j++){
        SI[i][j] = fac[j][i] / det_temp;
      }
    }
    for (i = 0; i < 3; i++){
      S[0][0][i][i] = SI[0][i];
      S[1][1][i][i] = SI[1][i];
      S[2][2][i][i] = SI[2][i];
      S[1][2][i][i] = SI[3][i] / 2.;
      S[2][0][i][i] = SI[4][i] / 2.;
      S[0][1][i][i] = SI[5][i] / 2.;
    }

    S[0][0][1][2] = SI[0][3] / 2.;
    S[1][1][1][2] = SI[1][3] / 2.;
    S[2][2][1][2] = SI[2][3] / 2.;
    S[1][2][1][2] = SI[3][3] / 4.;
    S[2][0][1][2] = SI[4][3] / 4.;
    S[0][1][1][2] = SI[5][3] / 4.;

    S[0][0][2][0] = SI[0][4] / 2.;
    S[1][1][2][0] = SI[1][4] / 2.;
    S[2][2][2][0] = SI[2][4] / 2.;
    S[1][2][2][0] = SI[3][4] / 4.;
    S[2][0][2][0] = SI[4][4] / 4.;
    S[0][1][2][0] = SI[5][4] / 4.;

    S[0][0][0][1] = SI[0][5] / 2.;
    S[1][1][0][1] = SI[1][5] / 2.;
    S[2][2][0][1] = SI[2][5] / 2.;
    S[1][2][0][1] = SI[3][5] / 4.;
    S[2][0][0][1] = SI[4][5] / 4.;
    S[0][1][0][1] = SI[5][5] / 4.;

    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {
              S[i][j][l][k] = S[i][j][k][l];
              S[j][i][k][l] = S[i][j][k][l];
              S[j][i][l][k] = S[i][j][k][l];
          }
        }
      }
    }
  } else {
    for (i=0; i<dimension; i++) {
      for (j=0; j<dimension; j++) {
        for (k=0; k<dimension; k++) {
          for (m=0; m<dimension; m++) {
            double isoS11, isoS12, isoS44;
            isoS11 = 1.0/(young);
            isoS12 = -xnu/(young);
            isoS44 = 2*(isoS11-isoS12);
            C[i][j][k][m] = mu * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+ll*DELTA(i,j)*DELTA(k,m);
            S[i][j][k][m] = isoS44/4 * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+isoS12*DELTA(i,j)*DELTA(k,m);
          }
        }
      }
    }
  }
  if (app->nphase>1) {
    int npa = app->nphase-1;
    int kpa;
      double CCDo[3][3][3][3], CCDI[6][6], SSDo[3][3][3][3], SSDI[6][6];
      double Cinc[npa][3][3][3][3];
      /* set Cijkl*/
      for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
          for (k=0; k<3; k++) {
            for (l=0; l<3; l++) {
              for ( kpa = 0; kpa < npa; kpa++) {
                Cinc[kpa][i][j][k][l] = 0.;
              }
              CCDo[i][j][k][l] = 0.;
              SSDo[i][j][k][l] = 0.;
            }
          }
        }
      }
      for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        CCDI[i][j] = 0.;
        SSDI[i][j] = 0.;
      }
      }
      if (aflag == 0) {
        for (i=0; i<6; i++) {
          for (j=0; j<6; j++) {
            CI[i][j] = 0.;
          }
        }
        // get property from material 2
        C11 = material2->C11;
        C12 = material2->C12;
        C44 = material2->C44;
        CI[0][0] = C11;
        CI[0][1] = C12;
        CI[3][3] = C44;
        CI[1][1] = CI[0][0];
        CI[2][2] = CI[0][0];
        CI[0][2] = CI[0][1];
        CI[1][2] = CI[0][1];
        CI[4][4] = CI[3][3];
        CI[5][5] = CI[3][3];
        for (i=0; i<6; i++) {
          for (j=0; j<6; j++) {
            if(fabs(CI[i][j]) > 1.E-5){ //JPC
              CI[j][i] = CI[i][j];
            }
          }
        }
        for (kpa = 0; kpa < npa; kpa++) {
          if (kpa == 0) {
            // no rotation
            for (i=0; i<6; i++) {
              for (j=0; j<6; j++) {
                  KCK[i][j] = CI[i][j];
              }
            }
            // BCC grid
            ed1[0] = prim[0][0]; ed1[1] =  prim[1][0]; ed1[2] = prim[2][0];
            ed2[0] = prim[0][1]; ed2[1] =  prim[1][1]; ed2[2] = prim[2][1];
            ed3[0] = prim[0][2]; ed3[1] =  prim[1][2]; ed3[2] = prim[2][2];
            rotate_stiffness(CI, ed1, ed2, ed3); 

            // left grain
             // ed1[0] = -1.0/sqrt(2.0); ed1[1] =  -1.0/sqrt(2.0); ed1[2] = 0.0;
             // ed2[0] = -1.0/sqrt(6.0); ed2[1] = 1.0/sqrt(6.0); ed2[2] = -2.0/sqrt(6.0);
             // ed3[0] = 1.0/sqrt(3.0); ed3[1] =  -1.0/sqrt(3.0); ed3[2] = -1.0/sqrt(3.0);
             // rotate_stiffness(CI, ed1, ed2, ed3);
             // GB matrix
             // KCK[0][0]=214.30; KCK[0][1]=98.35; KCK[0][2]=89.87; KCK[0][3]=1.97; KCK[0][4]=-0.53; KCK[0][5]=-3.45;
             // KCK[1][1]=214.82; KCK[1][2]=126.25; KCK[1][3]=-0.38; KCK[1][4]=-0.13; KCK[1][5]=-23.45;
             // KCK[2][2]=210.26; KCK[2][3]=0.39; KCK[2][4]=0.52; KCK[2][5]=-11.69;
             // KCK[3][3]=65.97; KCK[3][4]=-0.04; KCK[3][5]=-0.75;
             // KCK[4][4]=31.60; KCK[4][5]=0.0;
             // KCK[5][5]=44.81;
             // for (i=0; i<6; i++) {
             //   for (j=0; j<6; j++) {
             //       KCK[j][i] = KCK[i][j];
             //   }
             // }
             // for (i=0; i<6; i++) {
             //   for (j=0; j<6; j++) {
             //       KCK[i][j] = KCK[i][j]*1.E9;
             //   }
             // }
             // GB matrix
          }
          else {
            ed1[0] = prim[0][0]; ed1[1] =  prim[1][0]; ed1[2] = prim[2][0];
            ed2[0] = prim[0][1]; ed2[1] =  prim[1][1]; ed2[2] = prim[2][1];
            ed3[0] = prim[0][2]; ed3[1] =  prim[1][2]; ed3[2] = prim[2][2];
            rotate_stiffness(CI, ed1, ed2, ed3);
            // KCK[0][0]=214.30; KCK[0][1]=98.35; KCK[0][2]=89.87; KCK[0][3]=1.97; KCK[0][4]=-0.53; KCK[0][5]=-3.45;
            // KCK[1][1]=214.82; KCK[1][2]=126.25; KCK[1][3]=-0.38; KCK[1][4]=-0.13; KCK[1][5]=-23.45;
            // KCK[2][2]=210.26; KCK[2][3]=0.39; KCK[2][4]=0.52; KCK[2][5]=-11.69;
            // KCK[3][3]=67.97; KCK[3][4]=-0.04; KCK[3][5]=-0.75;
            // KCK[4][4]=31.60; KCK[4][5]=0.0;
            // KCK[5][5]=44.81;
            // for (i=0; i<6; i++) {
            //   for (j=0; j<6; j++) {
            //       KCK[j][i] = KCK[i][j];
            //   }
            // }
            // for (i=0; i<6; i++) {
            //   for (j=0; j<6; j++) {
            //       KCK[i][j] = KCK[i][j]*1E9;
            //   }
            // }

          }
          for (i = 0; i < 3; i++){
            Cinc[kpa][0][0][i][i] = KCK[0][i];
            Cinc[kpa][1][1][i][i] = KCK[1][i];
            Cinc[kpa][2][2][i][i] = KCK[2][i];
            Cinc[kpa][1][2][i][i] = KCK[3][i];
            Cinc[kpa][2][0][i][i] = KCK[4][i];
            Cinc[kpa][0][1][i][i] = KCK[5][i];
          }

          Cinc[kpa][0][0][1][2] = KCK[0][3];
          Cinc[kpa][1][1][1][2] = KCK[1][3];
          Cinc[kpa][2][2][1][2] = KCK[2][3];
          Cinc[kpa][1][2][1][2] = KCK[3][3];
          Cinc[kpa][2][0][1][2] = KCK[4][3];
          Cinc[kpa][0][1][1][2] = KCK[5][3];

          Cinc[kpa][0][0][2][0] = KCK[0][4];
          Cinc[kpa][1][1][2][0] = KCK[1][4];
          Cinc[kpa][2][2][2][0] = KCK[2][4];
          Cinc[kpa][1][2][2][0] = KCK[3][4];
          Cinc[kpa][2][0][2][0] = KCK[4][4];
          Cinc[kpa][0][1][2][0] = KCK[5][4];

          Cinc[kpa][0][0][0][1] = KCK[0][5];
          Cinc[kpa][1][1][0][1] = KCK[1][5];
          Cinc[kpa][2][2][0][1] = KCK[2][5];
          Cinc[kpa][1][2][0][1] = KCK[3][5];
          Cinc[kpa][2][0][0][1] = KCK[4][5];
          Cinc[kpa][0][1][0][1] = KCK[5][5];

          for (i=0; i<3; i++) {
            for (j=0; j<3; j++) {
              for (k=0; k<3; k++) {
                for (l=0; l<3; l++) {
                  if(fabs(Cinc[kpa][i][j][k][l]) > 1.E-5){
                    Cinc[kpa][i][j][l][k] = Cinc[kpa][i][j][k][l];
                    Cinc[kpa][j][i][k][l] = Cinc[kpa][i][j][k][l];
                    Cinc[kpa][j][i][l][k] = Cinc[kpa][i][j][k][l];
                  }
                }
              }
            }
          }

        }

        // for (i=0; i<6; i++) {
        //   for (j=0; j<6; j++) {
        //     KCK[i][j] = CI[i][j];
        //   }
        // }


      } else {
        for (kpa = 0; kpa < npa; kpa++) {
        double incC11, incC12, incC44, incmu, incll; 
        double sclC = 2.0;
        incC11 = sclC*C11;
        incC12 = sclC*C12;
        incC44 = sclC*C44;
       /* incC11 = 84.68E9;//249.01E9;//169.9E9;//246.5E9; //Mg inclusion
        incC12 = 56.6E9;//135.43E9;//122.6E9;//147.3E9;
        incC44 = 61.4E9;//18.1E9;//76.2E9;//127.4E9;*/
      /*  incC11 = 249.01E9;//84.68E9;//249.01E9;//169.9E9;//246.5E9; //Nb inclusion
        incC12 = 135.43E9;//56.6E9;//135.43E9;//122.6E9;//147.3E9;
        incC44 = 18.1E9;//61.4E9;//18.1E9;//76.2E9;//127.4E9; */
       /* incC11 = 168.4E9; //Cu inclusion
        incC12 = 121.4E9;
        incC44 = 75.4E9;*/

        // incmu = sclC*mu;//(incC11-incC12+3.0*incC44)/5.0;
        // incll = sclC*ll;//(incC11+2.0*incC12)/3.0 - 2.0/3.0*incmu;
        // (11,11,1)
        // incmu = 0.5*55.2E9;//38.53E9; // what is this? *** : inclusion mu and lambda (second phase) for the isotropic case
        // incll = 0.5*138.4E9;//108.5E9;
        // (12,12,1)
        //incmu = 37.87E9;
        //incll = 110.95E9
        // (2,2,5)
        // incmu = 35.61E9;
        // incll = 108.77E9;;
        // get properties from material 2
        incmu = material2->mu;
        incll = material2->ll;

        if(me == 0){
          printf("Inc mu %E, ll %E, incC12 %E\n", incmu, incll, incC12);
        }
        for (i=0; i<dimension; i++)
          for (j=0; j<dimension; j++)
            for (k=0; k<dimension; k++)
              for (m=0; m<dimension; m++) {
                 //Cinc[kpa][i][j][k][m] = incmu * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+incll*DELTA(i,j)*DELTA(k,m);
                Cinc[kpa][i][j][k][m] = incmu * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+incll*DELTA(i,j)*DELTA(k,m);
                //Cinc[kpa][i][j][k][m] = 23.5E9 * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+121.57E9*DELTA(i,j)*DELTA(k,m);
                //Cinc[kpa][i][j][k][m] = 92.03E9 * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+22.82E9*DELTA(i,j)*DELTA(k,m);
                //Cinc[kpa][i][j][k][m] = 42.44E9 * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+56.6E9*DELTA(i,j)*DELTA(k,m);
        }
      }
      }

      for (kpa = 0; kpa < npa; kpa++) {
        for (i=0; i<dimension; i++) {
          for (j=0; j<dimension; j++) {
            for (k=0; k<dimension; k++) {
              for (m=0; m<dimension; m++) {
                CCDo[i][j][k][m] = (Cinc[kpa][i][j][k][m] - C[i][j][k][m]) / 1.E9;
              }
            }
          }
        }
        for (i = 0; i < 3; i++){
          CCDI[0][i] = CCDo[0][0][i][i];
          CCDI[1][i] = CCDo[1][1][i][i];
          CCDI[2][i] = CCDo[2][2][i][i];
          CCDI[3][i] = CCDo[1][2][i][i];
          CCDI[4][i] = CCDo[2][0][i][i];
          CCDI[5][i] = CCDo[0][1][i][i];
        }

        CCDI[0][3] = CCDo[0][0][1][2];
        CCDI[1][3] = CCDo[1][1][1][2];
        CCDI[2][3] = CCDo[2][2][1][2];
        CCDI[3][3] = CCDo[1][2][1][2];
        CCDI[4][3] = CCDo[2][0][1][2];
        CCDI[5][3] = CCDo[0][1][1][2];

        CCDI[0][4] = CCDo[0][0][2][0];
        CCDI[1][4] = CCDo[1][1][2][0];
        CCDI[2][4] = CCDo[2][2][2][0];
        CCDI[3][4] = CCDo[1][2][2][0];
        CCDI[4][4] = CCDo[2][0][2][0];
        CCDI[5][4] = CCDo[0][1][2][0];

        CCDI[0][5] = CCDo[0][0][0][1];
        CCDI[1][5] = CCDo[1][1][0][1];
        CCDI[2][5] = CCDo[2][2][0][1];
        CCDI[3][5] = CCDo[1][2][0][1];
        CCDI[4][5] = CCDo[2][0][0][1];
        CCDI[5][5] = CCDo[0][1][0][1];

        k = 0;
        for (i=0; i<6; i++) {
          for (j = 0; j< 6; j++) {
            SSDI[i][j] = 0.;
            bb[i][j] = 0.;
            fac[i][j] = 0.;
            if(fabs(CCDI[i][j]) < 1.E-5){
              k = k+1;
            }
          }
        }
        // if(k == 36){
        //   if(rank == 0) {
        //     printf("CC in material %d equals that of material 1. \n", kdo+1);
        //   }
        //   exit(1);
        // }
        for (i=0; i<3; i++) {
          for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
              for (l=0; l<3; l++) {
                MM[kpa][i][j][k][l] = 0.;
              }
            }
          }
        }
        k = 6;
        det_temp = determinant(CCDI, k);
        // if(me == 0){
        //   printf("CCDI det= %E\n", det_temp);
        //   printf("CCDI 11 = %E; 44 = %E; 12 = %E;\n", CCDI[0][0],CCDI[3][3],CCDI[0][1]);
        // }

        for (q = 0;q < k; q++) {
          for (p = 0;p < k; p++) {
            m = 0;
            n = 0;
            for (i = 0;i < k; i++) {
              for (j = 0;j < k; j++){
                if (i != q && j != p){
                  bb[m][n] = CCDI[i][j];
                  if (n < (k - 2)) {
                    n++;
                  }
                  else {
                    n = 0;
                    m++;
                  }
                }
              }
            }
            fac[q][p] = pow(-1, q + p) * determinant(bb, k - 1);
          }
        }

        for (i = 0; i < k; i++){
          for (j = 0; j < k; j++){
            SSDI[i][j] = fac[j][i] / det_temp;
          }
        }
        // if(me == 0){
        //   printf("SSDI det= %E\n", det_temp);
        //   printf("SSDI 11 = %E; 44 = %E; 12 = %E;\n", SSDI[0][0],SSDI[3][3],SSDI[0][1]);
        // }

        for (i = 0; i < 3; i++){
          SSDo[0][0][i][i] = SSDI[0][i];
          SSDo[1][1][i][i] = SSDI[1][i];
          SSDo[2][2][i][i] = SSDI[2][i];
          SSDo[1][2][i][i] = SSDI[3][i];
          SSDo[2][0][i][i] = SSDI[4][i];
          SSDo[0][1][i][i] = SSDI[5][i];
        }

        SSDo[0][0][1][2] = SSDI[0][3];
        SSDo[1][1][1][2] = SSDI[1][3];
        SSDo[2][2][1][2] = SSDI[2][3];
        SSDo[1][2][1][2] = SSDI[3][3];
        SSDo[2][0][1][2] = SSDI[4][3];
        SSDo[0][1][1][2] = SSDI[5][3];

        SSDo[0][0][2][0] = SSDI[0][4];
        SSDo[1][1][2][0] = SSDI[1][4];
        SSDo[2][2][2][0] = SSDI[2][4];
        SSDo[1][2][2][0] = SSDI[3][4];
        SSDo[2][0][2][0] = SSDI[4][4];
        SSDo[0][1][2][0] = SSDI[5][4];

        SSDo[0][0][0][1] = SSDI[0][5];
        SSDo[1][1][0][1] = SSDI[1][5];
        SSDo[2][2][0][1] = SSDI[2][5];
        SSDo[1][2][0][1] = SSDI[3][5];
        SSDo[2][0][0][1] = SSDI[4][5];
        SSDo[0][1][0][1] = SSDI[5][5];

        for (i=0; i<3; i++) {
          for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
              for (l=0; l<3; l++) { // change units from 1/GPa to 1/Pa
                SSDo[i][j][k][l] = SSDo[i][j][k][l] / 1.E9;
              }
            }
          }
        }


        for (i=0; i<3; i++) {
          for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
              for (l=0; l<3; l++) {
                if(fabs(SSDo[k][l][i][j] - SSDo[i][j][k][l]) > 1.E-5){
                  if(me == 0){
                    printf("Fail 1st sym check SSDo1 = %lf, SSDo2 = %lf \n", SSDo[k][l][i][j], SSDo[i][j][k][l]);
                  }
                }
              }
            }
          }
        }

        for (i=0; i<3; i++) {
          for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
              for (l=0; l<3; l++) {
                if(fabs(SSDo[i][j][k][l]) > 1.E-5){
                  SSDo[i][j][l][k] = SSDo[i][j][k][l];
                  SSDo[j][i][k][l] = SSDo[i][j][k][l];
                  SSDo[j][i][l][k] = SSDo[i][j][k][l];
                }
              }
            }
          }
        }

        for (i=0; i<3; i++) {
          for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
              for (l=0; l<3; l++) {
                MM[kpa][i][j][k][l] = 0.;
                for (m=0; m<3; m++) {
                  for (n=0; n<3; n++) {
                    for (p=0; p<3; p++) {
                      for (q=0; q<3; q++) {
                        MM[kpa][i][j][k][l] += -C[i][j][m][n]*SSDo[m][n][p][q]*C[p][q][k][l]/mu;
                      }
                    }
                  }
                }
                MM[kpa][i][j][k][l] -= C[i][j][k][l]/mu; // normalized by mu[0]
              }
            }
          }
        }
        if(me == 0){
          printf("Mmatrix done M11 %E; M44 %E; M12 %E\n", MM[kpa][0][0][0][0],MM[kpa][1][2][1][2],MM[kpa][1][1][0][0]);
        }
      } //kpa

  } //if-nphase

}

/* ----------------------------------------------------------------------
   Calculate stiffness for CTB
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::stiffness_CTB()
{
  double C11, C12, C44, mu,ll,young,xnu;
  double CI[6][6], SI[6][6], bb[6][6], fac[6][6], ed1[3], ed2[3], ed3[3], det_temp;
  int i, j, k, m, p, q, n,l;
  C44 = material->C44;
  C12 = material->C12;
  C11 = material->C11;
  mu  = material->mu;
  ll = material->ll;
  young  = material->young;
  xnu = material->xnu;
  int aflag = material->aflag;

  /* set Cijkl*/
  if (aflag == 0) {
    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        CI[i][j] = 0.;
        KCK_CTB[i][j] = 0.;
      }
    }
    CI[0][0] = C11;
    CI[0][1] = C12;
    CI[3][3] = C44;
    CI[1][1] = CI[0][0];
    CI[2][2] = CI[0][0];
    CI[0][2] = CI[0][1];
    CI[1][2] = CI[0][1];
    CI[4][4] = CI[3][3];
    CI[5][5] = CI[3][3];
    for (i=0; i<6; i++) {
       for (j=0; j<6; j++) {
         if(fabs(CI[i][j]) > 1.E-5){ //JPC
           CI[j][i] = CI[i][j];
         }
       }
     }
     // no rotation
     for (i=0; i<6; i++) {
       for (j=0; j<6; j++) {
           KCK_CTB[i][j] = CI[i][j];
       }
     }

     // BCC grid
     /*ed3[0] = -1.0/sqrt(2.0); ed3[1] =  1.0/sqrt(2.0); ed3[2] = 0.0;
     ed2[0] = 1.0/sqrt(6.0); ed2[1] = 1.0/sqrt(6.0); ed2[2] = -2.0/sqrt(6.0);
     ed1[0] = 1.0/sqrt(3.0); ed1[1] =  1.0/sqrt(3.0); ed1[2] = 1.0/sqrt(3.0);
     rotate_stiffness(CI, ed1, ed2, ed3); */
    // CTB Grid Primitive vectors (ed matrix needs to have determinant of 1)
    // ed1[0] = -1.0/cbrt(6);ed1[1] = 0.0/cbrt(6);ed1[2] = -1.0/cbrt(6);
    // ed2[0] = -1.0/cbrt(6);ed2[1] = 1.0/cbrt(6);ed2[2] = 0.0/cbrt(6);
    // ed3[0] = 1.0/cbrt(6);ed3[1] = 4.0/cbrt(6);ed3[2] = -1.0/cbrt(6);

    ed1[0] = 1.0/cbrt(6);ed1[1] = 4.0/cbrt(6);ed1[2] = -1.0/cbrt(6);
    ed2[0] = -1.0/cbrt(6);ed2[1] = 1.0/cbrt(6);ed2[2] = 0.0/cbrt(6);
    ed3[0] = 0.0/cbrt(6);ed3[1] = 1.0/cbrt(6);ed3[2] = 1.0/cbrt(6);


     rotate_stiffness_CTB(CI, ed1, ed2, ed3);
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {
            C_CTB[i][j][k][l] = 0.;
            S_CTB[i][j][k][l] = 0.;
          }
        }
      }
    }

    for (i = 0; i < 3; i++){
      C_CTB[0][0][i][i] = KCK_CTB[0][i];
      C_CTB[1][1][i][i] = KCK_CTB[1][i];
      C_CTB[2][2][i][i] = KCK_CTB[2][i];
      C_CTB[1][2][i][i] = KCK_CTB[3][i];
      C_CTB[2][0][i][i] = KCK_CTB[4][i];
      C_CTB[0][1][i][i] = KCK_CTB[5][i];
    }

    C_CTB[0][0][1][2] = KCK_CTB[0][3];
    C_CTB[1][1][1][2] = KCK_CTB[1][3];
    C_CTB[2][2][1][2] = KCK_CTB[2][3];
    C_CTB[1][2][1][2] = KCK_CTB[3][3];
    C_CTB[2][0][1][2] = KCK_CTB[4][3];
    C_CTB[0][1][1][2] = KCK_CTB[5][3];

    C_CTB[0][0][2][0] = KCK_CTB[0][4];
    C_CTB[1][1][2][0] = KCK_CTB[1][4];
    C_CTB[2][2][2][0] = KCK_CTB[2][4];
    C_CTB[1][2][2][0] = KCK_CTB[3][4];
    C_CTB[2][0][2][0] = KCK_CTB[4][4];
    C_CTB[0][1][2][0] = KCK_CTB[5][4];

    C_CTB[0][0][0][1] = KCK_CTB[0][5];
    C_CTB[1][1][0][1] = KCK_CTB[1][5];
    C_CTB[2][2][0][1] = KCK_CTB[2][5];
    C_CTB[1][2][0][1] = KCK_CTB[3][5];
    C_CTB[2][0][0][1] = KCK_CTB[4][5];
    C_CTB[0][1][0][1] = KCK_CTB[5][5];

    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {
            if(fabs(C_CTB[i][j][k][l]) > 1.E-5){
              C_CTB[i][j][l][k] = C_CTB[i][j][k][l];
              C_CTB[j][i][k][l] = C_CTB[i][j][k][l];
              C_CTB[j][i][l][k] = C_CTB[i][j][k][l];
            }
          }
        }
      }
    }
    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        SI[i][j] = 0.;
      }
    }
    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        bb[i][j] = 0.;
        fac[i][j] = 0.;
      }
    }
    k = 6;

    det_temp = determinant(KCK_CTB, k);
    if(me == 0){
      printf("KCK_CTB det= %E\n", det_temp);
    }
    for (q = 0;q < k; q++) {
      for (p = 0;p < k; p++) {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++) {
          for (j = 0;j < k; j++){
            if (i != q && j != p){
              bb[m][n] = KCK_CTB[i][j];
              if (n < (k - 2)) {
                n++;
              }
              else {
                n = 0;
                m++;
              }
            }
          }
        }
        fac[q][p] = pow(-1, q + p) * determinant(bb, k - 1);
      }
    }

    for (i = 0; i < k; i++){
      for (j = 0; j < k; j++){
        SI[i][j] = fac[j][i] / det_temp;
        // if(me == 0)printf("SI[%d][%d] = %E\n",i,j,SI[i][j]);
      }
    }
    for (i = 0; i < 3; i++){
      S_CTB[0][0][i][i] = SI[0][i];
      S_CTB[1][1][i][i] = SI[1][i];
      S_CTB[2][2][i][i] = SI[2][i];
      S_CTB[1][2][i][i] = SI[3][i] / 2.;
      S_CTB[2][0][i][i] = SI[4][i] / 2.;
      S_CTB[0][1][i][i] = SI[5][i] / 2.;
    }

    S_CTB[0][0][1][2] = SI[0][3] / 2.;
    S_CTB[1][1][1][2] = SI[1][3] / 2.;
    S_CTB[2][2][1][2] = SI[2][3] / 2.;
    S_CTB[1][2][1][2] = SI[3][3] / 4.;
    S_CTB[2][0][1][2] = SI[4][3] / 4.;
    S_CTB[0][1][1][2] = SI[5][3] / 4.;

    S_CTB[0][0][2][0] = SI[0][4] / 2.;
    S_CTB[1][1][2][0] = SI[1][4] / 2.;
    S_CTB[2][2][2][0] = SI[2][4] / 2.;
    S_CTB[1][2][2][0] = SI[3][4] / 4.;
    S_CTB[2][0][2][0] = SI[4][4] / 4.;
    S_CTB[0][1][2][0] = SI[5][4] / 4.;

    S_CTB[0][0][0][1] = SI[0][5] / 2.;
    S_CTB[1][1][0][1] = SI[1][5] / 2.;
    S_CTB[2][2][0][1] = SI[2][5] / 2.;
    S_CTB[1][2][0][1] = SI[3][5] / 4.;
    S_CTB[2][0][0][1] = SI[4][5] / 4.;
    S_CTB[0][1][0][1] = SI[5][5] / 4.;

    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {
              S_CTB[i][j][l][k] = S_CTB[i][j][k][l];
              S_CTB[j][i][k][l] = S_CTB[i][j][k][l];
              S_CTB[j][i][l][k] = S_CTB[i][j][k][l];
              // if(me == 0)printf("S[%d][%d][%d][%d] = %E\n",i,j,k,l,S[i][j][k][l]);
          }
        }
      }
    }
  } else {
    for (i=0; i<dimension; i++) {
      for (j=0; j<dimension; j++) {
        for (k=0; k<dimension; k++) {
          for (m=0; m<dimension; m++) {
            double isoS11, isoS12, isoS44;
            isoS11 = 1.0/(young);
            isoS12 = -xnu/(young);
            isoS44 = 2*(isoS11-isoS12);
            C_CTB[i][j][k][m] = mu * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+ll*DELTA(i,j)*DELTA(k,m);
            S_CTB[i][j][k][m] = isoS44/4 * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+isoS12*DELTA(i,j)*DELTA(k,m);
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
Calculate Green's function
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::greens_function()
{
  // double Gaux[dimension][dimension];
  // double detG=0;
  // double fk[dimension];
  // double fk2, fk4;

  // for(int k1=0;k1<local_x;k1++){
  //   for(int k2=0;k2<local_y;k2++){
  //     for(int k3=0;k3<local_z;k3++){
  // 	int nfreq = k3+(k2)*local_z+(k1)*local_z*local_y;
  // 	fk[0] = fx[nfreq];
  // 	fk[1] = fy[nfreq];
  // 	fk[2] = fz[nfreq];
  // 	fk2 = fk[0]*fk[0]+fk[1]*fk[1]+fk[2]*fk[2];
  // 	fk4 = fk2*fk2;
  // 	if(fk2>0){
  // 	  for(int i=0; i<dimension; i++) {
  // 	    for(int k=0; k<dimension; k++) {
  // 	      Gaux[i][k] = 0.0;
  // 	      for(int j=0; j<dimension; j++) {
  // 		for(int l=0; l<dimension; l++) {
  // 		  Gaux[i][k] += C[k][j][i][l] * fk[l] * fk[j];
  // 		}
  // 	      }
  // 	    }
  // 	  }
  // 	}

  // 	// Now we have to compute the inverse
  // 	// This would be much better using Eigen

  // 	for(int i = 0; i < dimension; i++)
  // 	  detG += (Gaux[0][i] * (Gaux[1][(i+1)%dimension] * Gaux[2][(i+2)%dimension] - Gaux[1][(i+2)%dimension] * Gaux[2][(i+1)%dimension]));
  // 	for(int i = 0; i < 3; i++){
  // 	  for(int j = 0; j < 3; j++){
  //           G[k1][k2][k3][i][j] = ((Gaux[(j+1)%dimension][(i+1)%dimension] * Gaux[(j+2)%dimension][(i+2)%dimension])
  // 				   - (Gaux[(j+1)%dimension][(i+2)%dimension] * Gaux[(j+2)%dimension][(i+1)%dimension]))/ detG;
  // 	  }
  // 	}
  //     }
  //   }
  // }
}

/* ----------------------------------------------------------------------
generate BB matrix
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::Bmatrix()
{
  int i, j, k, l, m, n, u, v, k1, k2, k3, ka, kb, nv, nb, nfreq;
  int is, js, ks;
  double fkr;
  //double C[dimension][dimension][dimension][dimension];
  //double Crot[dimension][dimension][dimension][dimension];
  double A[dimension][dimension][dimension][dimension];
  // double B[slip_systems][slip_systems][local_x][local_y][local_z];
  double G[dimension][dimension], G_inv[dimension][dimension];
  //double Grot[dimension][dimension];
  double fk[dimension];
  double xn_temp[slip_systems][dimension], xb_temp[slip_systems][dimension], eps_temp[slip_systems][dimension][dimension];
  double xnu, fk2, fk4, fka, fkb, ll, mu, young, G_inv_det;
  double C44 = material->C44;
  double C12 = material->C12;
  double C11 = material->C11;
  double dslip = material->dslip;
  // int flag_aniso = materials->flag_aniso;
  double F[slip_systems][dimension][dimension];
  double D[slip_systems][dimension][dimension];
  double E[dimension][dimension][dimension][dimension];

  int ND = dimension;

  int aflag = material->aflag;
  int tflag = app->tflag;
  mu = material->mu;
  // material->ll = ll;
  young = material->young;
  // material->xnu = xnu;
  // mu = C44-(2.0*C44+C12-C11)/5.0;
  //ll = C12-(2.0*C44+C12-C11)/5.0;
  // young = mu*(3*ll+2*mu)/(ll+mu);
  ll = material->ll;
  xnu = material->xnu; //young/2.0/mu-1.0;

  if(me == 0){
    printf("Bmatrix mu %E, ll %E, young %E, nu %E\n", mu, ll, young, xnu);
  }

  for (int i=0; i<dimension; i++) {
    for (int j=0; j<dimension; j++) {
      for (int k=0; k<dimension; k++) {
        for (int m=0; m<dimension; m++) {
          A[i][j][k][m] = 0.0;
        }
      }
    }
  }

  for (i=0; i<slip_systems; i++) {
    for (j=0; j<dimension; j++) {
      xn_temp[i][j]=xn[i][j];
      xb_temp[i][j]=xb[i][j];
    }
  }

  /*set eps */
  // DO WE NEED THESE TWO LOOPS???
  for (k=0; k<slip_systems;k++){
    for (i=0; i<dimension; i++) {
      for (j=0; j<dimension; j++) {
          eps_temp[k][i][j] = xb_temp[k][i]*xn_temp[k][j]/dslip;
          eps[k][i][j] = eps_temp[k][i][j];
        }
      }
    }

  /* set A, Green function and B matrix*/
  for(k1=0;k1<local_x;k1++){
    for(k2=0;k2<local_y;k2++){
      for(k3=0;k3<local_z;k3++){
        nfreq = k3+(k2)*local_z+(k1)*local_z*local_y;
        fk[0] = fx[nfreq];
        fk[1] = fy[nfreq];
        fk[2] = fz[nfreq];
        fk2 = fk[0]*fk[0]+fk[1]*fk[1]+fk[2]*fk[2];
        fk4 = fk2*fk2;
        for (i=0; i<ND; i++) {
          for (k=0; k<ND; k++) {
            G[i][k] = 0.0;
          }
        }
        if(fk2>0){

          if (aflag == 0) {
            for (i=0; i<ND; i++) {
              for (k=0; k<ND; k++) {
                G_inv[i][k] = 0.0;
                for (j=0; j<ND; j++) {
                  for (l=0; l<ND; l++) {
                    if(tflag == 1 && ((k3 < int(local_z/8.0)) || k3 > int(local_z*7.0/8.0))){
                      G_inv[i][k] = G_inv[i][k] + C_CTB[i][j][k][l] * fk[j] * fk[l];
                    }
                    else{
                    G_inv[i][k] = G_inv[i][k] + C[i][j][k][l] * fk[j] * fk[l];
                  }
                }
              }
            }
            }

            G_inv_det = 0.;

            for(i = 0; i < 3; i++){  // determinant of G_inv
              G_inv_det = G_inv_det + (G_inv[0][i] * (G_inv[1][(i+1)%3] * G_inv[2][(i+2)%3] - G_inv[1][(i+2)%3] * G_inv[2][(i+1)%3]));
            }

            // calculate G from G_inv

            for(i = 0; i < 3; i++){
              for(j = 0; j < 3; j++) {
                G[j][i] = ((G_inv[(j+1)%3][(i+1)%3] * G_inv[(j+2)%3][(i+2)%3]) - (G_inv[(j+1)%3][(i+2)%3] * G_inv[(j+2)%3][(i+1)%3])) / G_inv_det;
              }
            }
          } else {

            for(k = 0; k < 3; k++){
              for(i = 0; i < 3; i++) {
                G[k][i] = (2.0 * DELTA(i,k)/fk2-1.0/(1.0-xnu)*fk[i]*fk[k]/fk4)/(2.0*mu);
              }
            }
          }

          for (m=0; m<dimension; m++) {
            for (n=0; n<dimension; n++) {
              for (u=0; u<dimension; u++) {
                for (v=0; v<dimension; v++) {
                  nb = nfreq + m*local_x*local_y*local_z + n*local_x*local_y*local_z*dimension + u*local_x*local_y*local_z*dimension*dimension + v*local_x*local_y*local_z*dimension*dimension*dimension;
                  AA[nb] = 0.0;
                  A[m][n][u][v] = 0.0;

                  for(i=0; i<dimension; i++) {
                    for(j=0; j<dimension; j++) {
                      for(k=0; k<dimension; k++) {
                        for(l=0; l<dimension; l++) {
                          if(tflag == 1 && ((k3 < int(local_z/8.0)) || k3 > int(local_z*7.0/8.0))){
                          A[m][n][u][v] = A[m][n][u][v] - C_CTB[k][l][u][v]*C_CTB[i][j][m][n]*G[k][i]*fk[j]*fk[l]/mu ;
                          }
                          else{
                          A[m][n][u][v] = A[m][n][u][v] - C[k][l][u][v]*C[i][j][m][n]*G[k][i]*fk[j]*fk[l]/mu ;
                        }
                      }
                    }
                  }
                  }
                  if(tflag == 1 && ((k3 < int(local_z/8.0)) || k3 > int(local_z*7.0/8.0))){
                    A[m][n][u][v] = A[m][n][u][v]+C_CTB[m][n][u][v]/mu;
                  }
                  else{
                    A[m][n][u][v] = A[m][n][u][v]+C[m][n][u][v]/mu;
                    AA[nb] = A[m][n][u][v];
                  }
                  /*printf("A %d %d %d %d %f %f %f %f\n",m,n,u,v,fk[0],fk[1],fk[2],A[m][n][u][v]);*/
                }
              }
            }
          }
        } /*if fk2 */
        for(ka=0;ka<slip_systems;ka++){
          for(kb=0;kb<slip_systems;kb++){
            nb = nfreq +(ka)*local_x*local_y*local_z+(kb)*local_x*local_y*local_z*slip_systems;
            BB[nb] = 0.0;
            for (m=0; m<dimension; m++){
              for (n=0; n<dimension; n++) {
                for (u=0; u<dimension; u++) {
                  for (v=0; v<dimension; v++) {
                    BB[nb] += A[m][n][u][v]*eps[ka][m][n]*eps[kb][u][v];
                  }
                }
              }
            }
            /*printf("%lf %lf %lf %lf \n", fx[nfreq], fy[nfreq], fz[nfreq], BB[nb]);*/
          } /*ka*/
        }/* kb*/

        for (ka=0; ka<slip_systems; ka++)
        for (i=0; i<dimension; i++)
        for (j=0; j<dimension; j++){
          nb = nfreq + (ka)*local_x*local_y*local_z + i*local_x*local_y*local_z*slip_systems + j*local_x*local_y*local_z*slip_systems*dimension;
          FF[nb] = 0.0;
          for(k=0; k<dimension; k++)
          for (l=0; l<dimension; l++){
            for (m=0; m<dimension; m++)
            for (n=0; n<dimension; n++){
              if(tflag == 1 && ((k3 < int(local_z/8.0)) || k3 > int(local_z*7.0/8.0))){
              FF[nb] += C_CTB[k][l][m][n]*G[j][k]*eps[ka][m][n]*fk[i]*fk[l];
              }
              else{
              FF[nb] += C[k][l][m][n]*G[j][k]*eps[ka][m][n]*fk[i]*fk[l];
            }
          }
        }
        }

        for (ka=0; ka<slip_systems; ka++) {
          for (u=0; u<dimension; u++) {
            for (v=0; v<dimension; v++){
              nb = nfreq + (ka)*local_x*local_y*local_z + u*local_x*local_y*local_z*slip_systems + v*local_x*local_y*local_z*slip_systems*dimension;
              DD[nb] = 0.0;
              for (m=0; m<dimension; m++) {
                for (n=0; n<dimension; n++) {
                  DD[nb] += A[m][n][u][v]*eps[ka][m][n];
                }
              }
            } //j
          } //i
        }


        if (app->nphase>1) {
          for (i=0; i<dimension; i++)
          for (j=0; j<dimension; j++){
            for (m=0; m<dimension; m++)
            for (n=0; n<dimension; n++){
              E[i][j][m][n]= 0.0;
              for(k=0; k<dimension; k++)
              for (l=0; l<dimension; l++){
              if(tflag == 1 && ((k3 < int(local_z/8.0)) || k3 > int(local_z*7.0/8.0))){
                E[i][j][m][n] = E[i][j][m][n] + C_CTB[k][l][m][n]*G[j][k]*fk[i]*fk[l];
              }
              else{
                E[i][j][m][n] = E[i][j][m][n] + C[k][l][m][n]*G[j][k]*fk[i]*fk[l];
              }
              }
              nb = nfreq + (i)*local_x*local_y*local_z + j*local_x*local_y*local_z*dimension + m*local_x*local_y*local_z*dimension*dimension + n*local_x*local_y*local_z*dimension*dimension*dimension;
              EE[nb] = E[i][j][m][n];
            }
          }
        }


      }	/*k1*/
    }	/*k2*/
  }	/*k3*/
  if(me == 0){
    printf("A B F D E matrix done\n");
  }


  return;
}

/* ----------------------------------------------------------------------
generate FF matrix
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::Fmatrix()
{
  int i, j, k, l, m, n, u, v, k1, k2, k3, ka, nv, nb, nfreq;
  int is, js, ks;
  double fkr;
  //double C[dimension][dimension][dimension][dimension];
  double F[slip_systems][dimension][dimension];
  double D[slip_systems][dimension][dimension];
  double G[dimension][dimension];
  double fk[dimension];
  double xnu, mu, ll, young, fk2, fk4, fka,fkb;
  double A[dimension][dimension][dimension][dimension];
  double E[dimension][dimension][dimension][dimension];

  double C44 = material->C44;
  double C12 = material->C12;
  double C11 = material->C11;
  double dslip = material->dslip;

  mu = material->mu;
  young = material->young;
  ll = C12-(2.0*C44+C12-C11)/5.0;
  xnu = young/2.0/mu-1.0;



  if(me == 0){
    printf("Fmatrix mu %lf, ll %lf, young %lf, nu %lf\n", mu, ll, young, xnu);
  }

  for (int i=0; i<dimension; i++) {
    for (int j=0; j<dimension; j++) {
      G[i][j] = 0.0;
      for (int k=0; k<dimension; k++) {
        for (int m=0; m<dimension; m++) {
          A[i][j][k][m] = 0.0;
        }
      }
    }
  }

  /* set Green function and F matrix*/

  for(k1=0;k1<local_x;k1++)
  for(k2=0;k2<local_y;k2++)
  for(k3=0;k3<local_z;k3++){
    nfreq = k3+(k2)*local_z+(k1)*local_z*local_y;
    fk[0] = fx[nfreq];
    fk[1] = fy[nfreq];
    fk[2] = fz[nfreq];
    fk2 = fk[0]*fk[0]+fk[1]*fk[1]+fk[2]*fk[2];
    fk4 = fk2*fk2;

    if(fk2>0){
      for(m=0; m<dimension; m++)
      for(n=0; n<dimension; n++)
      for(u=0; u<dimension; u++)
      for(v=0; v<dimension; v++){
        A[m][n][u][v] = 0.0;

        for(i=0; i<dimension; i++)
        for(j=0; j<dimension; j++)
        for (k=0; k<dimension; k++){
          G[k][i] = (2.0 * DELTA(i,k)/fk2-1.0/(1.0-xnu)*fk[i]*fk[k]/fk4)/(2.0*mu);
          for(l=0; l<dimension; l++){
            A[m][n][u][v] = A[m][n][u][v] - C[k][l][u][v]*C[i][j][m][n]*G[k][i]*fk[j]*fk[l] ;
          }
        }
        A[m][n][u][v] = A[m][n][u][v]+C[m][n][u][v];
      }
    } /*if fk2 */

    for (ka=0; ka<slip_systems; ka++)
    for (i=0; i<dimension; i++)
    for (j=0; j<dimension; j++){

      F[ka][i][j]= 0.0;
      D[ka][i][j]=0.0;

      for(k=0; k<dimension; k++)
      for (l=0; l<dimension; l++){
        for (m=0; m<dimension; m++)
        for (n=0; n<dimension; n++){
          F[ka][i][j] = F[ka][i][j] + C[k][l][m][n]*G[j][k]*eps[ka][m][n]*fk[i]*fk[l] ;
        }
        D[ka][i][j]=D[ka][i][j]+A[i][j][k][l]*eps[ka][k][l];
      }

      nb = nfreq + (ka)*local_x*local_y*local_z + i*local_x*local_y*local_z*slip_systems + j*local_x*local_y*local_z*slip_systems*dimension;
      FF[nb] = F[ka][i][j];

      nb = nfreq + (ka)*local_x*local_y*local_z + i*local_x*local_y*local_z*slip_systems + j*local_x*local_y*local_z*slip_systems*dimension;
      DD[nb] = D[ka][i][j];
    }

    for (i=0; i<dimension; i++)
    for (j=0; j<dimension; j++){
      for (m=0; m<dimension; m++)
      for (n=0; n<dimension; n++){
        E[i][j][m][n]= 0.0;
        for(k=0; k<dimension; k++)
        for (l=0; l<dimension; l++){
          E[i][j][m][n] = E[i][j][m][n] + C[k][l][m][n]*G[j][k]*fk[i]*fk[l];
        }
        nb = nfreq + (i)*local_x*local_y*local_z + j*local_x*local_y*local_z*dimension + m*local_x*local_y*local_z*dimension*dimension + n*local_x*local_y*local_z*dimension*dimension*dimension;
        EE[nb] = E[i][j][m][n];
      }
    }

  }/*k1,k2,k3*/
  return;
}
/* ----------------------------------------------------------------------
resolve shear stress
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::resolSS_Schmid()
{
  int is, i, j, k, l;
  int NS = slip_systems;
  int ND = dimension;
  double dslip = material->dslip;

  for (is=0;is<NS;is++){
    tau[is] = 0.0;
    for(i=0;i<ND;i++){
      for(j=0;j<ND;j++){
        tau[is] = tau[is]+(sigma_rot[i][j]*xn[is][j]*xb[is][i])/dslip;
        //tau[is] = tau[is]+(sigma[i][j]*xn[is][j]*xb[is][i])/dslip;
        /*printf("%d %d %lf %lf %lf\n", i,j,sigma[i][j],xn[is][j], xb[is][i]);*/
        /*
        if(me==0){
        if (logfile) fprintf(logfile, "is = %d i = %d j = %d sigma[i][j] = %lf xn[is][j] = %lf xb[is][i] = %lf\n", \
        is, i, j, sigma_rot[i][j], xn[is][j], xb[is][i]);
        if (screen) fprintf(screen, "is = %d i = %d j = %d sigma[i][j] = %lf xn[is][j] = %lf xb[is][i] = %lf\n", \
        is, i, j, sigma_rot[i][j], xn[is][j], xb[is][i]);
      }
      */

    }
  }
  if(me == 0){
    if (logfile) fprintf(logfile,"Schmid: Resolved shear stresses tau[%d] = %lf\n",is, tau[is]);
    if (screen) fprintf(screen,"Schmid: Resolved shear stresses tau[%d] = %lf\n",is, tau[is]);

  }
}
return;
}
/* ----------------------------------------------------------------------
resolve shear stress : MRSSP with angle
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::resolSS_non_Schmid()
{
  int is, i, j, k, l;
  int NS = slip_systems;
  int ND = dimension;
  double dslip = material->dslip;

  /* The non-glide streses are accounted for 'tau' calculation. (R. Groger et al., Acta Mater. 56 (2008)) */
  double S_tot, theta_n;
  double xn_mrss[3][3];  // maximum resolved shear stress plane which has angle of 'angle_to_110' to {110} plane
  double xnp[3][3]; // plane which has 60^0 to mrss plane
  double xnxb[3][3]; // xn_mrss x xb
  double xnpxb[3][3]; // xnp x xb

  double angle_to_110 = app->angle_to_110;
  double w1= app->w1 ;
  double w2= app->w2;
  double w3= app->w3;

  theta_n = 0.0;
  for (is=0;is<NS;is++){
    //    theta_n = acos(xn[is][2])*180./M_PI;

    /* mrss plane with the angle to 110 plane */
    xn_mrss[is][0]= -sin((theta_n+angle_to_110)*M_PI/180.);
    xn_mrss[is][1]= 0.0;
    xn_mrss[is][2]= cos((theta_n+angle_to_110)*M_PI/180.);

    /* plane which has 60^0 to {110} planes */
    xnp[is][0]= -sin((theta_n+angle_to_110+60)*M_PI/180.); // -sin(60)
    xnp[is][1]= 0.0;
    xnp[is][2]= cos((theta_n+angle_to_110+60)*M_PI/180.); // cos(60)

    theta_n+=120. ;
  }

  for (is=0;is<NS;is++){
    tau[is] = 0.0;
    /* n0 X b (cross products) for S2 term */
    xnxb[is][0] = xn_mrss[is][1]*xb[is][2]-xn_mrss[is][2]*xb[is][1];
    xnxb[is][1] = -xn_mrss[is][0]*xb[is][2]+xn_mrss[is][2]*xb[is][0];
    xnxb[is][2] = xn_mrss[is][0]*xb[is][1]-xn_mrss[is][1]*xb[is][0];

    /* n0' X b (cross product) for S3 term */
    xnpxb[is][0] = xnp[is][1]*xb[is][2]-xnp[is][2]*xb[is][1];
    xnpxb[is][1] = -xnp[is][0]*xb[is][2]+xnp[is][2]*xb[is][0];
    xnpxb[is][2] = xnp[is][0]*xb[is][1]-xnp[is][1]*xb[is][0];

    for(i=0;i<ND;i++){
      for(j=0;j<ND;j++){
        // sum the projection tensors for non glide stresses
        S_tot = xn_mrss[is][j]*xb[is][i] + w1*xnp[is][j]*xb[is][i] + w2*xn_mrss[is][j]*xnxb[is][i] + w3*xnp[is][j]*xnpxb[is][i];
        tau[is] = tau[is]+(sigma_rot[i][j]*S_tot)/dslip;

        /*printf("%d %d %lf %lf %lf\n", i,j,sigma[i][j],xn[is][j], xb[is][i]);*/
      }
    }
    if(me == 0){

      if (logfile) fprintf(logfile, "angle %lf, w1 %lf,w2 %lf,w3 %lf \n",angle_to_110,w1,w2,w3);
      if (screen) fprintf(screen, "angle %lf, w1 %lf,w2 %lf,w3 %lf \n",angle_to_110,w1,w2,w3);
      if (logfile) fprintf(logfile,"nonSchmid: Resolved shear stresses tau[%d] = %lf\n",is, tau[is]);
      if (screen) fprintf(screen,"nonSchmid: Resolved shear stresses tau[%d] = %lf\n",is, tau[is]);
      if (logfile) fprintf(logfile, "n X b = %lf %lf %lf\n",xnxb[is][0],xnxb[is][1],xnxb[is][2]);
      if (logfile) fprintf(logfile, "n' X b = %lf %lf %lf\n",xnpxb[is][0],xnpxb[is][1],xnpxb[is][2]);

    }

  }

  return;
}

/* ----------------------------------------------------------------------
MPI send function
------------------------------------------------------------------------- */

void FFTW_Slab_Conduct::fft_send(void *buf, int count, MPI_Datatype datatype, int dest,
  int tag, MPI_Comm comm)
  {
    MPI_Send(buf,count,datatype,dest,tag,comm);
  }

  /* ----------------------------------------------------------------------
  MPI receive function
  ------------------------------------------------------------------------- */

  void FFTW_Slab_Conduct::fft_recv(void *buf, int count, MPI_Datatype datatype, int dest,
    int tag, MPI_Comm comm, MPI_Status *status)
    {
      MPI_Recv(buf,count,datatype,dest,tag,comm,status);
    }
    /* ----------------------------------------------------------------------
    grad(xi) and line character angle theta
    ------------------------------------------------------------------------- */
    void FFTW_Slab_Conduct::gradient()
    {

      int isa, i, j, k, index0, index, na, nr, nl;
      int nxl, nxr, nyl, nyr, nzl, nzr, dny, dnz;
      int ig, ib;
      int rank, np;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &np);

      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NS = slip_systems;
      int nmpi = NS*N2*N3;
  double xi_r[nmpi], xi_l[nmpi], gradx_r[nmpi], gradx_l[nmpi], mpi_r[nmpi], mpi_l[nmpi]; //order parameters in the skin region
      double d1,d2,d3;
      d1 = d2 = d3 = 1.0;
  double tang[3], magg, magt;
  double ngrad[3];
      MPI_Status status;

      //initialize
      for (i=0; i<nmpi; i++) {
        xi_r[i] = 0.0;
        xi_l[i] = 0.0;
        mpi_r[i] = 0.0;
        mpi_l[i] = 0.0;
    gradx_r[i] = 0.0;
    gradx_l[i] = 0.0;
      }
      for (i=0; i<nmpi*lN1; i++) {
        gradx[i] = 0.0;
        grady[i] = 0.0;
        gradz[i] = 0.0;
        theta[i] = 0.0;
    nhn[i] = 0.0;
    // picked[i] = 0.0;
  }
  for(isa=0;isa<NS;isa++){
    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
	for(k=0;k<N3;k++){
	  index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
    // if (xo[index]==0.0 && picked[index]!= -1.0) {
    if (xo[index]<=0.0) { // LOCATION 1 to reset picked value based on barrier, <= -> == 
      picked[index] = 0.0;
    }

  }
      }
      tang[0] = 0.0;
      tang[1] = 0.0;
      tang[2] = 0.0;

      // We send the skin region to the neighboring domains
      //  It needs to work also when there is just 1 proc
      i = 0;
      for(isa=0;isa<NS;isa++){
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          nr = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
          index0 = j*N3 + k + isa*N2*N3;
          xi_r[index0] = xi[0][nr];
          if(np == 1)
          mpi_r[index0] = xi_r[index0];
        }
      }
      i = lN1-1;
      for(isa=0;isa<NS;isa++){
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          nl = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
          index0 = j*N3 + k + isa*N2*N3;
          xi_l[index0] = xi[0][nl];
          if(np == 1)
          mpi_l[index0] = xi_l[index0];
        }
      }

      if(np > 1){

        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == np-1) {
          MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
          MPI_Send(&xi_r, nmpi, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
        }
        else if(rank == 0){
          MPI_Send(&xi_r, nmpi, MPI_DOUBLE, np-1, 1, MPI_COMM_WORLD);
          MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
        }
        else{
          MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
          MPI_Send(&xi_r, nmpi, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) {
          MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, np-1, 0, MPI_COMM_WORLD, &status);
          MPI_Send(&xi_l, nmpi, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        }
        else if(rank == np-1){
          MPI_Send(&xi_l, nmpi, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
          MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
        }
        else{
          MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
          MPI_Send(&xi_l, nmpi, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
      }

      // Once we have the skin we can calculate the gradient

      for(isa=0;isa<NS;isa++){

        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          index0 = j*N3 + k + isa*N2*N3;
          na = 2 * index; // = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
          //gradx
          nxr = na + 2*N2*N3; // = 2*((i+1)*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
          nxl = na - 2*N2*N3; // = 2*((i-1)*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
          if (i==0) {
            gradx[index] = (xi[0][nxr]-mpi_l[index0])/2;
          }
          else if (i==lN1-1){
            gradx[index] = (mpi_r[index0]-xi[0][nxl])/2;
          }
          else {
            gradx[index] = (xi[0][nxr]-xi[0][nxl])/2;
          }

          //grady
          nyr = na + 2*N3; // = 2*(i*N2*N3 + (j+1)*N3 + k + isa*lN1*N2*N3);
          nyl = na - 2*N3; // = 2*(i*N2*N3 + (j-1)*N3 + k + isa*lN1*N2*N3);
          dny = 2*(N2-1)*N3;
          if (j==0) {
            grady[index] = (xi[0][nyr]-xi[0][na+dny])/2;
          }
          else if (j==N2-1){
            grady[index] = (xi[0][na-dny]-xi[0][nyl])/2;
          }
          else {
            grady[index] = (xi[0][nyr]-xi[0][nyl])/2;
          }

          //gradz
          nzr = na + 2; // = 2*(i*N2*N3 + j*N3 + (k+1) + isa*lN1*N2*N3);
          nzl = na - 2; // = 2*(i*N2*N3 + j*N3 + (k-1) + isa*lN1*N2*N3);
          dnz = 2*(N3-1);
          if (k==0) {
            gradz[index] = (xi[0][nzr]-xi[0][na+dnz])/2;
          }
          else if (k==N3-1){
            gradz[index] = (xi[0][na-dnz]-xi[0][nzl])/2;
          }
          else {
            gradz[index] = (xi[0][nzr]-xi[0][nzl])/2;
          }
        }//ijk
      }//isa

    // Hessian calculation

    i = 0;
    for(isa=0;isa<NS;isa++){
      for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
    nr = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
    index0 = j*N3 + k + isa*N2*N3;
    gradx_r[index0] = gradx[nr];
    if(np == 1)
      mpi_r[index0] = gradx_r[index0];
        }
    }
    i = lN1-1;
    for(isa=0;isa<NS;isa++){
      for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
    nl = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
    index0 = j*N3 + k + isa*N2*N3;
    gradx_l[index0] = gradx[nl];
    if(np == 1)
      mpi_l[index0] = gradx_l[index0];
        }
    }
    if(np > 1){

      MPI_Barrier(MPI_COMM_WORLD);

      if (rank == np-1) {
        MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&gradx_r, nmpi, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
      }
      else if(rank == 0){
        MPI_Send(&gradx_r, nmpi, MPI_DOUBLE, np-1, 1, MPI_COMM_WORLD);
        MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
      }
      else{
        MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&gradx_r, nmpi, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      if (rank == 0) {
        MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, np-1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&gradx_l, nmpi, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
      }
      else if(rank == np-1){
        MPI_Send(&gradx_l, nmpi, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
      }
      else{
        MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&gradx_l, nmpi, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    for(isa=0;isa<NS;isa++){
      for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      index = 6*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
      index0 = j*N3 + k + isa*N2*N3;
      na = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
      //dx
      nxr = na + N2*N3; // = (i+1)*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
      nxl = na - N2*N3; // = (i-1)*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
      if (i==0) {
        hess[index] = (gradx[nxr]-mpi_l[index0])/2.0;
      }
      else if (i==lN1-1){
        hess[index] = (mpi_r[index0]-gradx[nxl])/2.0;
      }
      else {
        hess[index] = (gradx[nxr]-gradx[nxl])/2.0;
      }
      //dy
      nyr = na + N3; // = i*N2*N3 + (j+1)*N3 + k + isa*lN1*N2*N3;
      nyl = na - N3; // = i*N2*N3 + (j-1)*N3 + k + isa*lN1*N2*N3;
      dny = (N2-1)*N3;
      if (j==0) {
        hess[index+3] = (gradx[nyr]-gradx[na+dny])/2.0;
        hess[index+1] = (grady[nyr]-grady[na+dny])/2.0;
        hess[index+5] = (gradz[nyr]-gradz[na+dny])/2.0;
      }
      else if (j==N2-1){
        hess[index+3] = (gradx[na-dny]-gradx[nyl])/2.0;
        hess[index+1] = (grady[na-dny]-grady[nyl])/2.0;
        hess[index+5] = (gradz[na-dny]-gradz[nyl])/2.0;
      }
      else {
        hess[index+3] = (gradx[nyr]-gradx[nyl])/2.0;
        hess[index+1] = (grady[nyr]-grady[nyl])/2.0;
        hess[index+5] = (gradz[nyr]-gradz[nyl])/2.0;
      }
      //dz
      nzr = na + 1; // = i*N2*N3 + j*N3 + (k+1) + isa*lN1*N2*N3;
      nzl = na - 1; // = i*N2*N3 + j*N3 + (k-1) + isa*lN1*N2*N3;
      dnz = N3-1;
      if (k==0) {
        hess[index+4] = (gradx[nzr]-gradx[na+dnz])/2.0;
        hess[index+2] = (gradz[nzr]-gradz[na+dnz])/2.0;
      }
      else if (k==N3-1){
        hess[index+4] = (gradx[na-dnz]-gradx[nzl])/2.0;
        hess[index+2] = (gradz[na-dnz]-gradz[nzl])/2.0;
      }
      else {
        hess[index+4] = (gradx[nzr]-gradx[nzl])/2.0;
        hess[index+2] = (gradz[nzr]-gradz[nzl])/2.0;
      }
    } // ijk
    } // isa

      // transformation of the gradient from the local to global grid
      double tempG[3];
      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          tempG[0] = gradx[index]; tempG[1] = grady[index]; tempG[2] = gradz[index];
          gradx[index] = tempG[0]*four[0][0] + tempG[1]*four[1][0] + tempG[2]*four[2][0];
          grady[index] = tempG[0]*four[0][1] + tempG[1]*four[1][1] + tempG[2]*four[2][1];
          gradz[index] = tempG[0]*four[0][2] + tempG[1]*four[1][2] + tempG[2]*four[2][2];
        }
      }
    // transformation of the Hessian from the local to global grid
    int m, n, u;
    double tempH[3][3], Rh[3][3], RhR[3][3];

    for(isa=0;isa<NS;isa++){
      for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
      for(k=0;k<N3;k++){
        index = 6*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
        for(m=0;m<3;m++)
        for(n=0;n<3;n++){
          tempH[m][n] = 0.0;
          Rh[m][n] = 0.0;
          RhR[m][n] = 0.0;
        }
        tempH[0][0] = hess[index];
        tempH[1][1] = hess[index+1];
        tempH[2][2] = hess[index+2];
        tempH[0][1] = hess[index+3];
        tempH[1][0] = hess[index+3];
        tempH[0][2] = hess[index+4];
        tempH[2][0] = hess[index+4];
        tempH[1][2] = hess[index+5];
        tempH[2][1] = hess[index+5];

        for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++)
        for (u = 0; u < 3; u++){
          Rh[m][n] += tempH[m][u]*four[u][n];
        }
        for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++)
        for (u = 0; u < 3; u++){
          RhR[m][n] += four[u][m]*Rh[u][n];
        }
        hess[index] = RhR[0][0];
        hess[index+1] = RhR[1][1];
        hess[index+2] = RhR[2][2];
        hess[index+3] = RhR[0][1];
        hess[index+4] = RhR[0][2];
        hess[index+5] = RhR[2][1];
      }
    }
  //calculate line character
      double dir;
      double dotp;
      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          magg = sqrt(gradx[index]*gradx[index]+grady[index]*grady[index]+gradz[index]*gradz[index]);
          if (magg<0.05 || (lxs+i)!=0) {
            theta[index]= -1; //define
          }
          else {
            ngrad[0] = gradx[index]/magg;
            ngrad[1] = grady[index]/magg;
            ngrad[2] = gradz[index]/magg;
            tang[0] = ngrad[1]*xn[isa][2] - ngrad[2]*xn[isa][1];
            tang[1] = ngrad[2]*xn[isa][0] - ngrad[0]*xn[isa][2];
            tang[2] = ngrad[0]*xn[isa][1] - ngrad[1]*xn[isa][0];
            magt = sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]);
            tang[0] /= magt;
            tang[1] /= magt;
            tang[2] /= magt;
            dotp = tang[0]*xb[isa][0]+tang[1]*xb[isa][1]+tang[2]*xb[isa][2];
            dotp = (dotp>-1.0) ? dotp : -1.0;
            dotp = (dotp<1.0) ? dotp : 1.0;
            theta[index] = acos(dotp);
            dir = (xb[isa][1]*tang[2]-xb[isa][2]*tang[1])*xn[isa][0]+
            (xb[isa][2]*tang[0]-xb[isa][0]*tang[2])*xn[isa][1]+
            (xb[isa][0]*tang[1]-xb[isa][1]*tang[0])*xn[isa][2];
            if (dir<0.0) {
              theta[index]=acos(-1.0)-theta[index];
            }
            // dir[index] = (xb[isa][1]*tang[2]-xb[isa][2]*tang[1])*xn[isa][0]+
            //       (xb[isa][2]*tang[0]-xb[isa][0]*tang[2])*xn[isa][1]+
            //       (xb[isa][0]*tang[1]-xb[isa][1]*tang[0])*xn[isa][2];
            // if (dir[index]<0.0) {
            //     theta[index]=acos(-1)-theta[index];
            // }
          }
        }//ijk
      }//isa

  for(isa=0;isa<NS;isa++){
    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
    index0 = 6*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
    magg = sqrt(gradx[index]*gradx[index]+grady[index]*grady[index]+gradz[index]*gradz[index]);
    if (magg>=0.05) {
    ngrad[0] = gradx[index]/magg;
    ngrad[1] = grady[index]/magg;
    ngrad[2] = gradz[index]/magg;
    nhn[index] = hess[index0]*ngrad[0]*ngrad[0]    +  hess[index0+1]*ngrad[1]*ngrad[1]  +  hess[index0+2]*ngrad[2]*ngrad[2]
             + 2*hess[index0+3]*ngrad[0]*ngrad[1] + 2*hess[index0+4]*ngrad[0]*ngrad[2] + 2*hess[index0+5]*ngrad[1]*ngrad[2];
    }
  } //ijk
  } //isa


  double cosx, cosy, cosz;
  int signx, signy, signz;
  int indexx,indexy,indexz;
  int nv = 1;
  int iv, indexv;
  int vertex[nv][dimension];
  double forward_sum = 0.0;
  double backward_sum = 0.0;
  double forward_mean = 0.0;
  double backward_mean = 0.0;
  int vcount = 0;
  double nvertex[nv][dimension];
  double s2 = sqrt(2);
  double s3 = sqrt(3);
  double angle[7];
  double maxang = 0;
  int maxid = 0;

    vertex[0][0] = 0; vertex[0][1] = 0; vertex[0][2] = 1;

    nvertex[0][0] = 0; nvertex[0][1] = 0; nvertex[0][2] = 1;

  // vertex[0][0] = 0; vertex[0][1] = 0; vertex[0][2] = 1;
  // vertex[1][0] = 0; vertex[1][1] = 1; vertex[1][2] = 0;
  // vertex[2][0] = 0; vertex[2][1] = 1; vertex[2][2] = 1;
  // vertex[3][0] = 1; vertex[3][1] = 0; vertex[3][2] = 0;
  // vertex[4][0] = 1; vertex[4][1] = 0; vertex[4][2] = 1;
  // vertex[5][0] = 1; vertex[5][1] = 1; vertex[5][2] = 0;
  // vertex[6][0] = 1; vertex[6][1] = 1; vertex[6][2] = 1;

  // nvertex[0][0] = 0;    nvertex[0][1] = 0;    nvertex[0][2] = 1;
  // nvertex[1][0] = 0;    nvertex[1][1] = 1;    nvertex[1][2] = 0;
  // nvertex[2][0] = 0;    nvertex[2][1] = 1/s2; nvertex[2][2] = 1/s2;
  // nvertex[3][0] = 1;    nvertex[3][1] = 0;    nvertex[3][2] = 0;
  // nvertex[4][0] = 1/s2; nvertex[4][1] = 0;    nvertex[4][2] = 1/s2;
  // nvertex[5][0] = 1/s2; nvertex[5][1] = 1/s2; nvertex[5][2] = 0;
  // nvertex[6][0] = 1/s3; nvertex[6][1] = 1/s3; nvertex[6][2] = 1/s3;

  for(isa=0;isa<NS;isa++){
    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
	for(k=0;k<N3;k++){
	  index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
    if (theta[index]>=0.0 && xo[index]<=0.0) { //Second xo condition, picked value ...
    magg = sqrt(gradx[index]*gradx[index]+grady[index]*grady[index]+gradz[index]*gradz[index]);
    ngrad[0] = gradx[index]/magg;
    ngrad[1] = grady[index]/magg;
    ngrad[2] = gradz[index]/magg;
    cosx = ngrad[0]*prim[0][0] + ngrad[1]*prim[1][0] + ngrad[2]*prim[2][0];
    cosy = ngrad[0]*prim[0][1] + ngrad[1]*prim[1][1] + ngrad[2]*prim[2][1];
    cosz = ngrad[0]*prim[0][2] + ngrad[1]*prim[1][2] + ngrad[2]*prim[2][2];
    signx = (cosx>=0) ? 1 : -1;
    signy = (cosy>=0) ? 1 : -1;
    signz = (cosz>=0) ? 1 : -1;
/* find closest forward/backward grid pts
    for (iv = 0; iv < nv; iv++) {
      angle[iv] = ngrad[0]*signx*nvertex[iv][0] + ngrad[1]*signy*nvertex[iv][1] + ngrad[2]*signz*nvertex[iv][2];
      angle[iv] = (angle[iv]>-1.0) ? angle[iv] : -1.0;
      angle[iv] = (angle[iv]<1.0) ? angle[iv] : 1.0;
    }
    maxang = angle[0];
    for (iv = 1; iv < nv; iv++) {
      if (maxang<angle[iv]) {
        maxang = (angle[iv]);
        maxid = iv;
      }
    }
    indexv = (i+signx*vertex[maxid][0])*N2*N3 + (j+signy*vertex[maxid][1])*N3 + (k+signz*vertex[maxid][2]) + isa*lN1*N2*N3;
    forward_mean = nhn[indexv];

    indexv = (i-signx*vertex[maxid][0])*N2*N3 + (j-signy*vertex[maxid][1])*N3 + (k-signz*vertex[maxid][2]) + isa*lN1*N2*N3;
    backward_mean = nhn[indexv];
*/

    for (iv = 0; iv < nv; iv++) {
      indexx = i+signx*vertex[iv][0];
      if (indexx<0) indexx = lN1-1;
      if (indexx==lN1) indexx = 0;
      indexy = j+signy*vertex[iv][1];
      if (indexy<0) indexy = N2-1;
      if (indexy==N2) indexy = 0;
      indexz = k+signz*vertex[iv][2];
      if(indexz<0) indexz = N3-1;
      if(indexz==N3) indexz = 0;
      indexv = indexx*N2*N3 + indexy*N3 + indexz + isa*lN1*N2*N3;
      if (theta_d[indexv]>=0.0){
        forward_sum += nhn[indexv];
        vcount += 1;
      }
    }
    forward_mean = (vcount == 0) ? 0.0 : (forward_sum / vcount);
    forward_sum = 0.0;
    vcount = 0;

    for (iv = 0; iv < nv; iv++) {
      indexx = i-signx*vertex[iv][0];
      if (indexx<0) indexx = lN1-1;
      if (indexx==lN1) indexx = 0;
      indexy = j-signy*vertex[iv][1];
      if (indexy<0) indexy = N2-1;
      if (indexy==N2) indexy = 0;
      indexz = k-signz*vertex[iv][2];
      if(indexz<0) indexz = N3-1;
      if(indexz==N3) indexz = 0;
      indexv = indexx*N2*N3 + indexy*N3 + indexz + isa*lN1*N2*N3;
      if (theta_d[indexv]>=0.0){
        backward_sum += nhn[indexv];
        vcount += 1;
      }
    }
    backward_mean = (vcount == 0) ? 0.0 : (backward_sum / vcount);
    backward_sum = 0.0;
    vcount = 0;

    if (forward_mean < 0.0 && backward_mean > 0.0) { //picked1
      picked[index] = 1.0;
    }
// dislocation center passed->reset kMC
    //if (xo[index]==-1.0 && picked[index]==0.0) xo[index] = 0.0;
    

    }
  } //ijk
  } //isa
}

/* ----------------------------------------------------------------------
   grad(xi) and line character angle theta for FCC
------------------------------------------------------------------------- */
void FFTW_Slab_Conduct::gradient_fcc()
{

  int isa, i, j, k, index0, index, index1, na, nr, nl;
  int nxl, nxr, nyl, nyr, nzl, nzr, dny, dnz;
  int ig, ib;
  int rank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NS = slip_systems;
  int NP = num_planes;
  int nmpi = NP*N2*N3;
  double xi_r[nmpi], xi_l[nmpi], gradx_r[nmpi], gradx_l[nmpi], mpi_r[nmpi], mpi_l[nmpi]; //order parameters in the skin region
  double d1,d2,d3;
  d1 = d2 = d3 = 1.0;
  double tang[3], magg, magt;
  double ngrad[3];
  MPI_Status status;

  //initialize
  for (i=0; i<nmpi; i++) {
    xi_r[i] = 0.0;
    xi_l[i] = 0.0;
    mpi_r[i] = 0.0;
    mpi_l[i] = 0.0;
    gradx_r[i] = 0.0;
    gradx_l[i] = 0.0;
  }
  for (i=0; i<nmpi*lN1; i++) {
    gradx_d[i] = 0.0;
    grady_d[i] = 0.0;
    gradz_d[i] = 0.0;
    theta_d[i] = 0.0;
    nhn_d[i] = 0.0;
    // picked[i] = 0.0;
  }
  for(isa=0;isa<NP;isa++){
    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
	for(k=0;k<N3;k++){
	  index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
    index1 = i*N2*N3 + j*N3 + k + (3*isa)*lN1*N2*N3;
    // if (xo[index]==0.0 && picked[index]!= -1.0) {
    // if (xo[index1]<=0.0) { //changed from <= 0.0
    picked_d[index] = 0.0;
    // }

  }
  }
  tang[0] = 0.0;
  tang[1] = 0.0;
  tang[2] = 0.0;

  // We send the skin region to the neighboring domains
  //  It needs to work also when there is just 1 proc
  i = 0;
  for(isa=0;isa<NP;isa++){
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++){
	nr = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
	index0 = j*N3 + k + isa*N2*N3;
	xi_r[index0] = dreg[0][nr];
	if(np == 1)
	  mpi_r[index0] = xi_r[index0];
      }
  }
  i = lN1-1;
  for(isa=0;isa<NP;isa++){
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++){
	nl = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
	index0 = j*N3 + k + isa*N2*N3;
	xi_l[index0] = dreg[0][nl];
	if(np == 1)
	  mpi_l[index0] = xi_l[index0];
      }
  }

  if(np > 1){

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&xi_r, nmpi, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&xi_r, nmpi, MPI_DOUBLE, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&xi_r, nmpi, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
      MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, np-1, 0, MPI_COMM_WORLD, &status);
      MPI_Send(&xi_l, nmpi, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
    }
    else if(rank == np-1){
      MPI_Send(&xi_l, nmpi, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
      MPI_Send(&xi_l, nmpi, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // Once we have the skin we can calculate the gradient

  for(isa=0;isa<NP;isa++){

    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
	for(k=0;k<N3;k++){
	  index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
	  index0 = j*N3 + k + isa*N2*N3;
	  na = 2 * index; // = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
	  //gradx
	  nxr = na + 2*N2*N3; // = 2*((i+1)*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
	  nxl = na - 2*N2*N3; // = 2*((i-1)*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
	  if (i==0) {
	    gradx_d[index] = (dreg[0][nxr]-mpi_l[index0])/2;
	  }
	  else if (i==lN1-1){
	    gradx_d[index] = (mpi_r[index0]-dreg[0][nxl])/2;
	  }
	  else {
	    gradx_d[index] = (dreg[0][nxr]-dreg[0][nxl])/2;
	  }

	  //grady
	  nyr = na + 2*N3; // = 2*(i*N2*N3 + (j+1)*N3 + k + isa*lN1*N2*N3);
	  nyl = na - 2*N3; // = 2*(i*N2*N3 + (j-1)*N3 + k + isa*lN1*N2*N3);
	  dny = 2*(N2-1)*N3;
	  if (j==0) {
	    grady_d[index] = (dreg[0][nyr]-dreg[0][na+dny])/2;
	  }
	  else if (j==N2-1){
	    grady_d[index] = (dreg[0][na-dny]-dreg[0][nyl])/2;
	  }
	  else {
	    grady_d[index] = (dreg[0][nyr]-dreg[0][nyl])/2;
	  }

	  //gradz
	  nzr = na + 2; // = 2*(i*N2*N3 + j*N3 + (k+1) + isa*lN1*N2*N3);
	  nzl = na - 2; // = 2*(i*N2*N3 + j*N3 + (k-1) + isa*lN1*N2*N3);
	  dnz = 2*(N3-1);
	  if (k==0) {
	    gradz_d[index] = (dreg[0][nzr]-dreg[0][na+dnz])/2;
	  }
	  else if (k==N3-1){
	    gradz_d[index] = (dreg[0][na-dnz]-dreg[0][nzl])/2;
	  }
	  else {
	    gradz_d[index] = (dreg[0][nzr]-dreg[0][nzl])/2;
	  }
	}//ijk
  }//isa

  // Hessian calculation

  i = 0;
  for(isa=0;isa<NP;isa++){
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++){
  nr = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
  index0 = j*N3 + k + isa*N2*N3;
  gradx_r[index0] = gradx_d[nr];
  if(np == 1)
    mpi_r[index0] = gradx_r[index0];
      }
  }
  i = lN1-1;
  for(isa=0;isa<NP;isa++){
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++){
  nl = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
  index0 = j*N3 + k + isa*N2*N3;
  gradx_l[index0] = gradx_d[nl];
  if(np == 1)
    mpi_l[index0] = gradx_l[index0];
      }
  }
  if(np > 1){

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&gradx_r, nmpi, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&gradx_r, nmpi, MPI_DOUBLE, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&mpi_r, nmpi, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&gradx_r, nmpi, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
      MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, np-1, 0, MPI_COMM_WORLD, &status);
      MPI_Send(&gradx_l, nmpi, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
    }
    else if(rank == np-1){
      MPI_Send(&gradx_l, nmpi, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&mpi_l, nmpi, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
      MPI_Send(&gradx_l, nmpi, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for(isa=0;isa<NP;isa++){
    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
	for(k=0;k<N3;k++){
    index = 6*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
    index0 = j*N3 + k + isa*N2*N3;
	  na = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
	  //dx
	  nxr = na + N2*N3; // = (i+1)*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
	  nxl = na - N2*N3; // = (i-1)*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
	  if (i==0) {
	    hess_d[index] = (gradx_d[nxr]-mpi_l[index0])/2.0;
	  }
	  else if (i==lN1-1){
	    hess_d[index] = (mpi_r[index0]-gradx_d[nxl])/2.0;
	  }
	  else {
	    hess_d[index] = (gradx_d[nxr]-gradx_d[nxl])/2.0;
	  }
    //dy
    nyr = na + N3; // = i*N2*N3 + (j+1)*N3 + k + isa*lN1*N2*N3;
    nyl = na - N3; // = i*N2*N3 + (j-1)*N3 + k + isa*lN1*N2*N3;
    dny = (N2-1)*N3;
    if (j==0) {
	    hess_d[index+3] = (gradx_d[nyr]-gradx_d[na+dny])/2.0;
      hess_d[index+1] = (grady_d[nyr]-grady_d[na+dny])/2.0;
      hess_d[index+5] = (gradz_d[nyr]-gradz_d[na+dny])/2.0;
	  }
	  else if (j==N2-1){
	    hess_d[index+3] = (gradx_d[na-dny]-gradx_d[nyl])/2.0;
      hess_d[index+1] = (grady_d[na-dny]-grady_d[nyl])/2.0;
      hess_d[index+5] = (gradz_d[na-dny]-gradz_d[nyl])/2.0;
	  }
	  else {
	    hess_d[index+3] = (gradx_d[nyr]-gradx_d[nyl])/2.0;
      hess_d[index+1] = (grady_d[nyr]-grady_d[nyl])/2.0;
      hess_d[index+5] = (gradz_d[nyr]-gradz_d[nyl])/2.0;
	  }
    //dz
	  nzr = na + 1; // = i*N2*N3 + j*N3 + (k+1) + isa*lN1*N2*N3;
	  nzl = na - 1; // = i*N2*N3 + j*N3 + (k-1) + isa*lN1*N2*N3;
	  dnz = N3-1;
	  if (k==0) {
	    hess_d[index+4] = (gradx_d[nzr]-gradx_d[na+dnz])/2.0;
      hess_d[index+2] = (gradz_d[nzr]-gradz_d[na+dnz])/2.0;
	  }
	  else if (k==N3-1){
	    hess_d[index+4] = (gradx_d[na-dnz]-gradx_d[nzl])/2.0;
      hess_d[index+2] = (gradz_d[na-dnz]-gradz_d[nzl])/2.0;
	  }
	  else {
	    hess_d[index+4] = (gradx_d[nzr]-gradx_d[nzl])/2.0;
      hess_d[index+2] = (gradz_d[nzr]-gradz_d[nzl])/2.0;
	  }
  } // ijk
  } // isa

  // transformation of the gradient from the local to global grid
  double tempG[3];
  for(isa=0;isa<NP;isa++){
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
      tempG[0] = gradx_d[index]; tempG[1] = grady_d[index]; tempG[2] = gradz_d[index];
      gradx_d[index] = tempG[0]*four[0][0] + tempG[1]*four[1][0] + tempG[2]*four[2][0];
      grady_d[index] = tempG[0]*four[0][1] + tempG[1]*four[1][1] + tempG[2]*four[2][1];
      gradz_d[index] = tempG[0]*four[0][2] + tempG[1]*four[1][2] + tempG[2]*four[2][2];
    }
  }
  // transformation of the Hessian from the local to global grid
  int m, n, u;
  double tempH[3][3], Rh[3][3], RhR[3][3];

  for(isa=0;isa<NP;isa++){
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      index = 6*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
      for(m=0;m<3;m++)
      for(n=0;n<3;n++){
        tempH[m][n] = 0.0;
        Rh[m][n] = 0.0;
        RhR[m][n] = 0.0;
      }
      tempH[0][0] = hess_d[index];
      tempH[1][1] = hess_d[index+1];
      tempH[2][2] = hess_d[index+2];
      tempH[0][1] = hess_d[index+3];
      tempH[1][0] = hess_d[index+3];
      tempH[0][2] = hess_d[index+4];
      tempH[2][0] = hess_d[index+4];
      tempH[1][2] = hess_d[index+5];
      tempH[2][1] = hess_d[index+5];

      for (m = 0; m < 3; m++)
      for (n = 0; n < 3; n++)
      for (u = 0; u < 3; u++){
        Rh[m][n] += tempH[m][u]*four[u][n];
      }
      for (m = 0; m < 3; m++)
      for (n = 0; n < 3; n++)
      for (u = 0; u < 3; u++){
        RhR[m][n] += four[u][m]*Rh[u][n];
      }
      hess_d[index] = RhR[0][0];
      hess_d[index+1] = RhR[1][1];
      hess_d[index+2] = RhR[2][2];
      hess_d[index+3] = RhR[0][1];
      hess_d[index+4] = RhR[0][2];
      hess_d[index+5] = RhR[2][1];
    }
  }
  //calculate line character
  double dir;
  double dotp;
  for(isa=0;isa<NS;isa++){
    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
    magg = sqrt(gradx_d[index]*gradx_d[index]+grady_d[index]*grady_d[index]+gradz_d[index]*gradz_d[index]);
    if ((lxs+i)!=0) {
      theta_d[index]= -1; //define
    }
  //   else {
  //     ngrad[0] = gradx[index]/magg;
  //     ngrad[1] = grady[index]/magg;
  //     ngrad[2] = gradz[index]/magg;
  //     tang[0] = ngrad[1]*xn[isa][2] - ngrad[2]*xn[isa][1];
  //     tang[1] = ngrad[2]*xn[isa][0] - ngrad[0]*xn[isa][2];
  //     tang[2] = ngrad[0]*xn[isa][1] - ngrad[1]*xn[isa][0];
  //     magt = sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]);
  //     tang[0] /= magt;
  //     tang[1] /= magt;
  //     tang[2] /= magt;
  //     dotp = tang[0]*xb[isa][0]+tang[1]*xb[isa][1]+tang[2]*xb[isa][2];
  //     dotp = (dotp>-1.0) ? dotp : -1.0;
  //     dotp = (dotp<1.0) ? dotp : 1.0;
  //     theta[index] = acos(dotp);
  //     dir = (xb[isa][1]*tang[2]-xb[isa][2]*tang[1])*xn[isa][0]+
  //       (xb[isa][2]*tang[0]-xb[isa][0]*tang[2])*xn[isa][1]+
  //       (xb[isa][0]*tang[1]-xb[isa][1]*tang[0])*xn[isa][2];
  //     if (dir<0.0) {
  //       theta[index]=acos(-1.0)-theta[index];
  //     }
  //   }
  }//ijk
  }//isa

  for(isa=0;isa<NP;isa++){
    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
    index0 = 6*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
    magg = sqrt(gradx_d[index]*gradx_d[index]+grady_d[index]*grady_d[index]+gradz_d[index]*gradz_d[index]);
    if (magg>=0.05) {
    ngrad[0] = gradx_d[index]/magg;
    ngrad[1] = grady_d[index]/magg;
    ngrad[2] = gradz_d[index]/magg;
    nhn_d[index] = hess_d[index0]*ngrad[0]*ngrad[0]    +  hess_d[index0+1]*ngrad[1]*ngrad[1]  +  hess_d[index0+2]*ngrad[2]*ngrad[2]
             + 2*hess_d[index0+3]*ngrad[0]*ngrad[1] + 2*hess_d[index0+4]*ngrad[0]*ngrad[2] + 2*hess_d[index0+5]*ngrad[1]*ngrad[2];
    }
  } //ijk
  } //isa


  double cosx, cosy, cosz;
  int signx, signy, signz;
  int nv = 3;
  int iv, indexv;
  int indexx,indexy,indexz;
  int vertex[nv][dimension];
  double forward_sum = 0.0;
  double backward_sum = 0.0;
  double forward_mean = 0.0;
  double backward_mean = 0.0;
  int vcount = 0;
  double nvertex[nv][dimension];
  double s2 = sqrt(2);
  double s3 = sqrt(3);
  double s5 = sqrt(5);
  double angle[7];
  double maxang = 0;
  int maxid = 0;

  vertex[0][0] = 0; vertex[0][1] = 0; vertex[0][2] = 1;
  vertex[1][0] = 0; vertex[1][1] = -1; vertex[1][2] = 1;
  vertex[2][0] = 0; vertex[2][1] = -1; vertex[2][2] = 2;

  nvertex[0][0] = 0; nvertex[0][1] = 0; nvertex[0][2] = 1;
  nvertex[1][0] = 0; nvertex[1][1] = -1/s2; nvertex[1][2] = 1/s2;
  nvertex[2][0] = 0; nvertex[2][1] = -1/s5; nvertex[2][2] = 2/s5;

  // vertex[0][0] = 0; vertex[0][1] = 0; vertex[0][2] = 1;
  // vertex[1][0] = 0; vertex[1][1] = 1; vertex[1][2] = 0;
  // vertex[2][0] = 0; vertex[2][1] = 1; vertex[2][2] = 1;
  // vertex[3][0] = 1; vertex[3][1] = 0; vertex[3][2] = 0;
  // vertex[4][0] = 1; vertex[4][1] = 0; vertex[4][2] = 1;
  // vertex[5][0] = 1; vertex[5][1] = 1; vertex[5][2] = 0;
  // vertex[6][0] = 1; vertex[6][1] = 1; vertex[6][2] = 1;

  // nvertex[0][0] = 0;    nvertex[0][1] = 0;    nvertex[0][2] = 1;
  // nvertex[1][0] = 0;    nvertex[1][1] = 1;    nvertex[1][2] = 0;
  // nvertex[2][0] = 0;    nvertex[2][1] = 1/s2; nvertex[2][2] = 1/s2;
  // nvertex[3][0] = 1;    nvertex[3][1] = 0;    nvertex[3][2] = 0;
  // nvertex[4][0] = 1/s2; nvertex[4][1] = 0;    nvertex[4][2] = 1/s2;
  // nvertex[5][0] = 1/s2; nvertex[5][1] = 1/s2; nvertex[5][2] = 0;
  // nvertex[6][0] = 1/s3; nvertex[6][1] = 1/s3; nvertex[6][2] = 1/s3;

  for(isa=0;isa<NP;isa++){
    for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
	for(k=0;k<N3;k++){
	  index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
	  index1 = i*N2*N3 + j*N3 + k + (3*isa)*lN1*N2*N3;
    if (theta_d[index]>=0) {
    magg = sqrt(gradx_d[index]*gradx_d[index]+grady_d[index]*grady_d[index]+gradz_d[index]*gradz_d[index]);
    ngrad[0] = gradx_d[index]/magg;
    ngrad[1] = grady_d[index]/magg;
    ngrad[2] = gradz_d[index]/magg;
    cosx = ngrad[0]*prim[0][0] + ngrad[1]*prim[1][0] + ngrad[2]*prim[2][0];
    cosy = ngrad[0]*prim[0][1] + ngrad[1]*prim[1][1] + ngrad[2]*prim[2][1];
    cosz = ngrad[0]*prim[0][2] + ngrad[1]*prim[1][2] + ngrad[2]*prim[2][2];
    signx = (cosx>=0) ? 1 : -1;
    signy = (cosy>=0) ? 1 : -1;
    signz = (cosz>=0) ? 1 : -1;
/* find closest forward/backward grid pts
    for (iv = 0; iv < nv; iv++) {
      angle[iv] = ngrad[0]*signx*nvertex[iv][0] + ngrad[1]*signy*nvertex[iv][1] + ngrad[2]*signz*nvertex[iv][2];
      angle[iv] = (angle[iv]>-1.0) ? angle[iv] : -1.0;
      angle[iv] = (angle[iv]<1.0) ? angle[iv] : 1.0;
    }
    maxang = angle[0];
    for (iv = 1; iv < nv; iv++) {
      if (maxang<angle[iv]) {
        maxang = (angle[iv]);
        maxid = iv;
      }
    }
    indexv = (i+signx*vertex[maxid][0])*N2*N3 + (j+signy*vertex[maxid][1])*N3 + (k+signz*vertex[maxid][2]) + isa*lN1*N2*N3;
    forward_mean = nhn[indexv];

    indexv = (i-signx*vertex[maxid][0])*N2*N3 + (j-signy*vertex[maxid][1])*N3 + (k-signz*vertex[maxid][2]) + isa*lN1*N2*N3;
    backward_mean = nhn[indexv];
*/  
    for (iv = 0; iv < nv; iv++) {
      indexx = i+signx*vertex[iv][0];
      if (indexx<0) indexx = lN1-1;
      if (indexx==lN1) indexx = 0;
      indexy = j+signy*vertex[iv][1];
      if (indexy<0) indexy = N2-1;
      if (indexy==N2) indexy = 0;
      indexz = k+signz*vertex[iv][2];
      if(indexz<0) indexz = N3-1;
      if(indexz==N3) indexz = 0;
      indexv = indexx*N2*N3 + indexy*N3 + indexz + isa*lN1*N2*N3;
      if (theta_d[indexv]>=0.0){
        forward_sum += nhn_d[indexv];
        vcount += 1;
      }
    }
    forward_mean = (vcount == 0) ? 0.0 : (forward_sum / vcount);
    forward_sum = 0.0;
    vcount = 0;

    for (iv = 0; iv < nv; iv++) {
      indexx = i-signx*vertex[iv][0];
      if (indexx<0) indexx = lN1-1;
      if (indexx==lN1) indexx = 0;
      indexy = j-signy*vertex[iv][1];
      if (indexy<0) indexy = N2-1;
      if (indexy==N2) indexy = 0;
      indexz = k-signz*vertex[iv][2];
      if(indexz<0) indexz = N3-1;
      if(indexz==N3) indexz = 0;
      indexv = indexx*N2*N3 + indexy*N3 + indexz + isa*lN1*N2*N3;
      if (theta_d[indexv]>=0.0){
        backward_sum += nhn_d[indexv];
        vcount += 1;
      }
    }
    backward_mean = (vcount == 0) ? 0.0 : (backward_sum / vcount);
    backward_sum = 0.0;
    vcount = 0;

    if (forward_mean < 0.0 && backward_mean > 0.0) { //picked1
      picked_d[index] = 1.0;
    }
// dislocation center passed->reset kMC
    //if (xo[index]==-1.0 && picked[index]==0.0) xo[index] = 0.0;
    

    }
  } //ijk
  } //isa

}
/* -----------------------------------------------------------------------
   KMC events
   ---------------------------------------------------------------------*/
// void FFTW_Slab_Conduct::events(index)
// {
  /*
  if (kflag == 1) {
    if (picked[index]<=1) {
      xo[index] = 0.0;
      picked[index] = -1.0;
    }
    else picked[index] -= 1;
  } */

/*
  for(int i=0; i<app->NE; i++){
    app->event_cond[i] = 0;
    app->event_time[i] = -1.0;
  }
  app->kflag = 0;
  app->cur_time = app->t_kmc;
  */
// }
    /* -----------------------------------------------------------------------
    bcc core energy for perfect dislocations
    ---------------------------------------------------------------------*/
    void FFTW_Slab_Conduct::core_energy_bcc_perfect()
    {
      int i, j, k, isa, index, num;
      int tag;

      int ND = dimension;
      int N1 = nx;
      double size = static_cast<double>(N1);
      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NP = num_planes;
      int NS = slip_systems;
      double dslip = material->dslip;

      int NT = solve->max_iter;
      int NSI = app->stoptime;

      double usf = material->usf;
      double b = material->b;
      double mu = material->mu;
      //   if(me == 0){
      //     if (logfile) fprintf(logfile,"b = %e\n", b);
      //     if (screen) fprintf(screen,"b = %e, usf = %lf, mu = %lf\n", b,usf,mu);
      // }
      // double ll = material->ll;
      // double young = material->young;
      // double xnu = material->xnu;
      usf = usf/(mu*dslip*b);

      // Checking without core energy
      //An = 0;
      double cof;
      E_core = 0;

      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          if (theta[index]<0.0) {
            cof = 1.0;
          }
          else if( theta[index]<1.23 ) {
            cof= 1.1603*theta[index]*theta[index] - 2.0431*theta[index]+1;
          }

          else{
            cof=0.5473*theta[index]*theta[index] - 2.0035*theta[index]+1.8923;
          }
          E_core += usf*cof*(sin(M_PI*data_fftw[index].re)*sin(M_PI*data_fftw[index].re))/N1;   //To make it general fft->data_fftw hasz to be general
          dE_core[index] = usf*cof*M_PI*sin(2.0*M_PI*data_fftw[index].re);
        }/*ijk*/
      }/*isa*/
    }

    /* -----------------------------------------------------------------------
    core energy for perfect dislocations
    ---------------------------------------------------------------------*/
    void FFTW_Slab_Conduct::core_energy_perfect()
    {
      int i, j, k, isa, index, num, index0;
      int tag;

      int ND = dimension;
      int N1 = nx;
      double size = static_cast<double>(N1);
      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NP = num_planes;
      int NS = slip_systems;
      double dslip = material->dslip;

      int NT = solve->max_iter;
      int NSI = app->stoptime;

      double c0 = material->c0;
      double c1 = material->c1;
      double c2 = material->c2;
      double c3 = material->c3;
      double c4 = material->c4;
      double a1 = material->a1;
      double a3 = material->a3;
      double isf = material->isf;
      double usf = material->usf;
      double An = material->An;
      double Cn = material->Cn;
      double b = material->b;

      double mu = material->mu;
      double ll = material->ll;
      double young = material->young;
      double xnu = material->xnu;

      c0 = c0/(mu);
      c1 = c1/(mu);
      c2 = c2/(mu);
      c3 = c3/(mu);
      c4 = c4/(mu);
      a1 = a1/(mu);
      a3 = a3/(mu);
      isf = isf/(mu*dslip*b);
      An = An/(mu*dslip*b);
      Cn = (usf - (isf/2.0))/(mu*dslip*b);

      E_core = 0;
      // temp 2phase core energy
     double An1;
     //An1 = 158.0E-3/(mu*dslip*b);
     //An1 = 158.0E-3/40.42*38.53/(mu*dslip*b);//459.89E-3/(mu*dslip*b); //Mg // This value is from slide 8 in PFDD-BCCMg_Janel ***
     //An1 = An;
     An = 158.0E-3/(mu*dslip*b);//676.775E-3/(mu*dslip*b); //Nb // Its hard coded here... *** : hard coded can be changed :)
     An1 = 158.0E-3/(mu*dslip*b);
    //  An1 = 1.9*An; //(40.42/38.53)*An; //original
      //An1 = 620.4E-3/(mu*dslip*b); //Cu
      //An = 933.4E-3/(mu*dslip*b); //Ni
     //if(me == 0) printf("An %E, An1 %E\n", An, An1);
     

      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          E_core += An*(sin(M_PI*data_fftw[index].re)*sin(M_PI*data_fftw[index].re))/N1;   //To make it general fft->data_fftw hasz to be general
          dE_core[index] = An*M_PI*sin(2.0*M_PI*data_fftw[index].re);

          if (app->nphase>1) {
            index0 = i*N2*N3 + j*N3 + k;
            // if (phase[index0] == 1) {
              dE_core[index] = An1*M_PI*sin(2.0*M_PI*data_fftw[index].re);
            // }
          }

        }/*ijk*/
      }/*isa*/
    }
/* -----------------------------------------------------------------------
   fcc core energy for dislocation intersections
   ---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_test()
{
  int i, j, k, isa, index, index1, index2, index3, index4, index5, index6, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1,xi2,xi3,xi4,xi5,xi6;
  double c0 = material->c0;
  double c1 = material->c1;
  double c2 = material->c2;
  double c3 = material->c3;
  double c4 = material->c4;
  double a1 = material->a1;
  double a3 = material->a3;
  double isf = material->isf;
  double usf = material->usf;
  double An = material->An;
  double Cn = material->Cn;
  double b = material->b;

  double mu = material->mu;

  c0 = c0/(mu*dslip*b);
  c1 = c1/(mu*dslip*b);
  c2 = c2/(mu*dslip*b);
  c3 = c3/(mu*dslip*b);
  c4 = c4/(mu*dslip*b);
  a1 = a1/(mu*dslip*b);
  a3 = a3/(mu*dslip*b);
  isf = isf/(mu*dslip*b);
  An = An/(mu*dslip*b);

  for(isa=0;isa<NS;isa++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3;
          index4 = i*N2*N3 + j*N3 + k + 3*lN1*N2*N3;
          index5 = i*N2*N3 + j*N3 + k + 4*lN1*N2*N3;
          index6 = i*N2*N3 + j*N3 + k + 5*lN1*N2*N3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;
          xi4 = data_fftw[index4].re;
          xi5 = data_fftw[index5].re;
          xi6 = data_fftw[index6].re;

          fcore[index] = An*(1-cos(M_PI*(xi2-xi3+xi5-xi6))*cos(M_PI*(xi1-xi3+xi6-xi4))*cos(M_PI*(xi1-xi2-xi4+xi5)));

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = -An*M_PI*cos(M_PI*(xi2 - xi3 + xi5 - xi6))*sin(M_PI*(-2*xi1 + xi2 + xi3 + 2*xi4 - xi5 - xi6));
          /*partial derivative wrt phase field 2*/
          dE_core[index2] = -An*M_PI*cos(M_PI*(xi1 - xi3 - xi4 + xi6))*sin(M_PI*(xi1 - 2*xi2 + xi3 - xi4 + xi6));
          /*partial derivative wrt phase field 3*/
          dE_core[index3] = -An*M_PI*cos(M_PI*(xi1 - xi2 - xi4 + xi5))*sin(M_PI*(xi1 + xi2 - 2*xi3 - xi4 + xi5));
          /*partial derivative wrt phase field 4*/
          dE_core[index4] = An*M_PI*cos(M_PI*(xi2 - xi3 + xi5 - xi6))*sin(M_PI*(-2*xi1 + xi2 + xi3 + 2*xi4 - xi5 - xi6));
          /*partial derivative wrt phase field 5*/
          dE_core[index5] = -An*M_PI*cos(M_PI*(xi1 - xi3 - xi4 + xi6))*sin(M_PI*(-1*xi1 + xi3 + xi4 - 2*xi5 + xi6));
          /*partial derivative wrt phase field 6*/
          dE_core[index6] = An*M_PI*cos(M_PI*(xi1 - xi2 - xi4 + xi5))*sin(M_PI*(xi1 - xi2 - xi4 - xi5 + 2*xi6));

        }/*ijk*/
      }
    }
  }/*plane*/
}
    /* -----------------------------------------------------------------------
    core energy for variable USFE
    ---------------------------------------------------------------------*/
    void FFTW_Slab_Conduct::core_energy_mpea()
    {
  // int i, j, k, isa, index, num, core_index;
  // int tag;

  // int ND = dimension;
  // int N1 = nx;
  // double size = static_cast<double>(N1);
  // int lN1 = local_x;
  // int lxs = local_x_start;
  // int N2 = local_y;
  // int N3 = local_z;
  // int NP = num_planes;
  // int NS = slip_systems;
  // double dslip = material->dslip;

  // int NT = solve->max_iter;
  // int NSI = app->stoptime;

  // double * A_mpea = material->A_mpea;
  // double b = material->b;

  // double mu = material->mu;
  // double ll = material->ll;
  // double young = material->young;
  // double xnu = material->xnu;

  // E_core = 0;

  // for(isa=0;isa<NS;isa++){
  //   for(i=0;i<lN1;i++)
  //     for(j=0;j<N2;j++)
  //       for(k=0;k<N3;k++){
  //         index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
  //         core_index = (lxs+i)*N2*N3 + j*N3 + k + isa*N1*N2*N3;
  //         E_core += A_mpea[core_index]/(mu*dslip*b)*(sin(M_PI*data_fftw[index].re)*sin(M_PI*data_fftw[index].re))/N1;
  //         dE_core[index] = A_mpea[core_index]/(mu*dslip*b)*M_PI*sin(2.0*M_PI*data_fftw[index].re);
  //   }/*ijk*/
  // }/*isa*/

  int i, j, k, u, isa, index, index1, index2, index3, num, core_index;
      int tag;

      int ND = dimension;
      int N1 = nx;
      double size = static_cast<double>(N1);
      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NP = num_planes;
      int NS = slip_systems;
      double dslip = material->dslip;

      int NT = solve->max_iter;
      int NSI = app->stoptime;

  double * fcc7_mpea = material->fcc7_mpea;
  double c0,c1,c2,c3,c4,a1,a3;
      double b = material->b;

      double mu = material->mu;
      double ll = material->ll;
      double young = material->young;
      double xnu = material->xnu;

  double xi1,xi2,xi3;

      E_core = 0;


  for(isa=0;isa<NP;isa++){
        for(i=0;i<lN1;i++)
          for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
              index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          core_index = (lxs+i)*N2*N3 + j*N3 + k + isa*7*N1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + isa*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + isa*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + isa*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          c0 = fcc7_mpea[0*N1*N2*N3 + core_index]/(mu*dslip*b);
          c1 = fcc7_mpea[1*N1*N2*N3 + core_index]/(mu*dslip*b);
          c2 = fcc7_mpea[2*N1*N2*N3 + core_index]/(mu*dslip*b);
          c3 = fcc7_mpea[3*N1*N2*N3 + core_index]/(mu*dslip*b);
          c4 = fcc7_mpea[4*N1*N2*N3 + core_index]/(mu*dslip*b);
          a1 = fcc7_mpea[5*N1*N2*N3 + core_index]/(mu*dslip*b);
          a3 = fcc7_mpea[6*N1*N2*N3 + core_index]/(mu*dslip*b);

          fcore[index] = (c0 + c1*(cos(2.0*M_PI*(xi1-xi2)) + cos(2.0*M_PI*(xi2-xi3)) + cos(2.0*M_PI*(xi3-xi1))) + c2*(cos(2.0*M_PI*(2.0*xi1-xi2-xi3)) + cos(2.0*M_PI*(2.0*xi2-xi3-xi1)) + cos(2.0*M_PI*(2.0*xi3-xi1-xi2))) + c3*(cos(4.0*M_PI*(xi1-xi2)) + cos(4.0*M_PI*(xi2-xi3)) + cos(4.0*M_PI*(xi3-xi1))) + c4*(cos(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + cos(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + cos(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + cos(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + cos(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + cos(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1*(sin(2.0*M_PI*(xi1-xi2)) + sin(2.0*M_PI*(xi2-xi3)) + sin(2.0*M_PI*(xi3-xi1))) + a3*(sin(4.0*M_PI*(xi1-xi2)) + sin(4.0*M_PI*(xi2-xi3)) + sin(4.0*M_PI*(xi3-xi1))));

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = 2.0*M_PI*(c1*(sin(2.0*M_PI*(xi2-xi1)) + sin(2.0*M_PI*(xi3-xi1))) + c2*(2.0*sin(2.0*M_PI*(xi2+xi3-2.0*xi1)) + sin(2.0*M_PI*(2.0*xi2-xi3-xi1)) + sin(2.0*M_PI*(2.0*xi3-xi1-xi2))) + (2.0*c3)*(sin(4.0*M_PI*(xi2-xi1)) + sin(4.0*M_PI*(xi3-xi1))) + c4*(3.0*sin(2.0*M_PI*(xi2+2.0*xi3-3.0*xi1)) + 3.0*sin(2.0*M_PI*(2.0*xi2+xi3-3.0*xi1)) + 2.0*sin(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + sin(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + sin(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + 2.0*sin(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1*(cos(2.0*M_PI*(xi1-xi2)) - cos(2.0*M_PI*(xi3-xi1))) + (2.0*a3)*(cos(4.0*M_PI*(xi1-xi2)) - cos(4.0*M_PI*(xi3-xi1))));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = 2.0*M_PI*(c1*(sin(2.0*M_PI*(xi1-xi2)) + sin(2.0*M_PI*(xi3-xi2))) + c2*(sin(2.0*M_PI*(2.0*xi1-xi2-xi3)) + 2.0*sin(2.0*M_PI*(xi3+xi1-2.0*xi2)) + sin(2.0*M_PI*(2.0*xi3-xi1-xi2))) + (2.0*c3)*(sin(4.0*M_PI*(xi1-xi2)) + sin(4.0*M_PI*(xi3-xi2))) + c4*(sin(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + 2.0*sin(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + 3.0*sin(2.0*M_PI*(xi3+2.0*xi1-3.0*xi2)) + 3.0*sin(2.0*M_PI*(2.0*xi3+xi1-3.0*xi2)) + 2.0*sin(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + sin(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1*(cos(2.0*M_PI*(xi2-xi3)) - cos(2.0*M_PI*(xi1-xi2))) + (2.0*a3)*(cos(4.0*M_PI*(xi2-xi3)) - cos(4.0*M_PI*(xi1-xi2))));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = 2.0*M_PI*(c1*(sin(2.0*M_PI*(xi2-xi3)) + sin(2.0*M_PI*(xi1-xi3))) + c2*(sin(2.0*M_PI*(2.0*xi1-xi2-xi3)) + sin(2.0*M_PI*(2.0*xi2-xi3-xi1)) + 2.0*sin(2.0*M_PI*(xi1+xi2-2.0*xi3))) + (2.0*c3)*(sin(4.0*M_PI*(xi2-xi3)) + sin(4.0*M_PI*(xi1-xi3))) + c4*(2.0*sin(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + sin(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + sin(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + 2.0*sin(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + 3.0*sin(2.0*M_PI*(xi1+2.0*xi2-3.0*xi3)) + 3.0*sin(2.0*M_PI*(2.0*xi1+xi2-3.0*xi3))) + a1*(cos(2.0*M_PI*(xi3-xi1)) - cos(2.0*M_PI*(xi2-xi3))) + (2.0*a3)*(cos(4.0*M_PI*(xi3-xi1)) - cos(4.0*M_PI*(xi2-xi3))));

        }/*ijk*/
      }/*isa*/
    }

void FFTW_Slab_Conduct::core_energy_lup2D()
{
  int i, j, k, ia, ja, ka, isa, index, index1, index2, index3, plane, num, s;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  int Nx = material->Nx_lup - 1;
  int Ny = material->Ny_lup - 1;
  double mu = material->mu;
  double b = material->b;
  double xi1, xi2, xi3;
  double op_pfdd[3];
  double opd_pfdd[3];
  double opx,opy,opdx,opdy;
  double x1, x2, y1, y2, f11, f12, f21, f22;

  double dE_core_test = 0.0;

  double * sfe_dx = material->sfe_dx;
  double * sfe_dy = material->sfe_dy;

  for(plane=0;plane < NP; plane++) {
    for(ia=0;ia<lN1;ia++) {
      for(ja=0;ja<N2;ja++) {
        for(ka=0;ka<N3;ka++) {

          index1 = ia*N2*N3 + ja*N3 + ka + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = ia*N2*N3 + ja*N3 + ka + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = ia*N2*N3 + ja*N3 + ka + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          op_pfdd[0] = data_fftw[index1].re;
          op_pfdd[1] = data_fftw[index2].re;
          op_pfdd[2] = data_fftw[index3].re;

          /* conversion */

          opx = op_pfdd[0] - 0.5 * (op_pfdd[1] + op_pfdd[2]);
          opy = 0.5 * (op_pfdd[1] - op_pfdd[2]);
          if ((opx < 0.0) || (opx >= 1.0)) {
            opx = opx - floor(opx);
          }
          if ((opy < 0.0) || (opy >= 1.0)){
            opy = opy - floor(opy);
          }

          /* interpolation */

          i = (int)(opx*Nx);
          j = (int)(opy*Ny);

          x1 = (double)i/Nx;
          x2 = (double)(i+1)/Nx;
          y1 = (double)j/Ny;
          y2 = (double)(j+1)/Ny;

          //sfe_dx
          if ((i == 0) && (j == 0)){
            f11 = sfe_dx[(Ny-1)*Nx+(Nx-1)];
            f12 = sfe_dx[j*Nx+(Nx-1)];
            f21 = sfe_dx[(Ny-1)*Nx+i];
            f22 = sfe_dx[j*Nx+i];
          }
          else if (i == 0) {
            f11 = sfe_dx[(j-1)*Nx+(Nx-1)];
            f12 = sfe_dx[j*Nx+(Nx-1)];
            f21 = sfe_dx[(j+1)*Nx+i];
            f22 = sfe_dx[j*Nx+i];
          }
          else if (j == 0) {
            f11 = sfe_dx[(Ny-1)*Nx+(i-1)];
            f12 = sfe_dx[j*Nx+(i-1)];
            f21 = sfe_dx[(Ny-1)*Nx+i];
            f22 = sfe_dx[j*Nx+i];
          }
          else {
            f11 = sfe_dx[(j-1)*Nx+(i-1)];
            f12 = sfe_dx[j*Nx+(i-1)];
            f21 = sfe_dx[(j-1)*Nx+i];
            f22 = sfe_dx[j*Nx+i];
          }

          opdx = 1. / ((x2-x1)*(y2-y1)) * (
            f11 * (x2-opx) * (y2-opy) +
            f21 * (opx-x1) * (y2-opy) +
            f12 * (x2-opx) * (opy-y1) +
            f22 * (opx-x1) * (opy-y1));

          // if(opdx != 0 && opx > 0.95)
          // printf("TEST %1.10f %f %f %f %d %d %1.10f %1.10f\n",opdx,opx,x1,x2,i,j,sfe_dx[i][1].re,sfe_dx[i-1][1].re);
          //sfe_dy
          if ((i == 0) && (j == 0)){
            f11 = sfe_dy[(Ny-1)*Nx+(Nx-1)];
            f12 = sfe_dy[j*Nx+(Nx-1)];
            f21 = sfe_dy[(Ny-1)*Nx+i];
            f22 = sfe_dy[j*Nx+i];
          }
          else if (i == 0) {
            f11 = sfe_dy[(j-1)*Nx+(Nx-1)];
            f12 = sfe_dy[j*Nx+(Nx-1)];
            f21 = sfe_dy[(j+1)*Nx+i];
            f22 = sfe_dy[j*Nx+i];
          }
          else if (j == 0) {
            f11 = sfe_dy[(Ny-1)*Nx+(i-1)];
            f12 = sfe_dy[j*Nx+(i-1)];
            f21 = sfe_dy[(Ny-1)*Nx+i];
            f22 = sfe_dy[j*Nx+i];
          }
          else {
            f11 = sfe_dy[(j-1)*Nx+(i-1)];
            f12 = sfe_dy[j*Nx+(i-1)];
            f21 = sfe_dy[(j-1)*Nx+i];
            f22 = sfe_dy[j*Nx+i];
          }

          opdy = 1. / ((x2-x1)*(y2-y1)) * (
            f11 * (x2-opx) * (y2-opy) +
            f21 * (opx-x1) * (y2-opy) +
            f12 * (x2-opx) * (opy-y1) +
            f22 * (opx-x1) * (opy-y1));

            /* conversion */
            opd_pfdd[0] = opdx;
            opd_pfdd[1] = -0.5 * opdx + 0.5 * opdy;
            opd_pfdd[2] = -0.5 * opdx - 0.5 * opdy;
            for(s = 0; s < /*NS/NP*/3; s++){
              index = ia*N2*N3 + ja*N3 + ka + s*lN1*N2*N3 + plane*lN1*N2*N3*3;
              dE_core[index] = opd_pfdd[s];
              // if(/*(lxs+i) == 0 && */((ka <= int(N3/8.0)) || ka >= int(N3*7.0/8.0))){
              //   dE_core[index] = dE_core[index] * 5;
              // }
              // if((lxs+ia) == 0 && plane == 0) dE_core_test += abs(dE_core[index]);
              }
            }
          }
        }//ka
      }//plane
      // if(me == 0) printf("AVERAGE dE_core: %lf\n",dE_core_test/(3*N2*N3));
}
/* -----------------------------------------------------------------------
core energy 1D lookup table
---------------------------------------------------------------------*/

void FFTW_Slab_Conduct::core_energy_lup1D()
{
  int i, j, k, ia, ja, ka, isa, index, s , num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NS = slip_systems;
  int Nx = material->Nx_lup - 1;
  int Ny = material->Ny_lup - 1;
  double op_pfdd[3];
  double opd_pfdd[3];
  double opx,opy,opdx;
  double x1, x2, y1, y2, f11, f12, f21, f22;

  double * sfe_dx = material->sfe_dx;

  for(s=0;s < NS; s++) {
    for(ia=0;ia<lN1;ia++) {
      for(ja=0;ja<N2;ja++) {
        for(ka=0;ka<N3;ka++) {

          index = ia*N2*N3 + ja*N3 + ka + s*lN1*N2*N3;
          op_pfdd[s] = data_fftw[index].re;

          opx = op_pfdd[s];
          if ((opx < 0.0) || (opx > 1.0)){
            opx = opx - floor(opx);
          }


          /* interpolation */

          i = (int)(opx*Nx);
          j = (int)(opy*Ny);

          x1 = (double)i/Nx;
          x2 = (double)(i+1)/Nx;
          y1 = (double)j/Ny;
          y2 = (double)(j+1)/Ny;

          //sfe_dx
          if ((i == 0) && (j == 0)){
            f11 = sfe_dx[(Ny-1)*Nx+(Nx-1)];
            f12 = sfe_dx[j*Nx+(Nx-1)];
            f21 = sfe_dx[(Ny-1)*Nx+i];
            f22 = sfe_dx[j*Nx+i];
          }
          else if (i == 0) {
            f11 = sfe_dx[(j-1)*Nx+(Nx-1)];
            f12 = sfe_dx[j*Nx+(Nx-1)];
            f21 = sfe_dx[(j+1)*Nx+i];
            f22 = sfe_dx[j*Nx+i];
          }
          else if (j == 0) {
            f11 = sfe_dx[(Ny-1)*Nx+(i-1)];
            f12 = sfe_dx[j*Nx+(i-1)];
            f21 = sfe_dx[(Ny-1)*Nx+i];
            f22 = sfe_dx[j*Nx+i];
          }
          else {
            f11 = sfe_dx[(j-1)*Nx+(i-1)];
            f12 = sfe_dx[j*Nx+(i-1)];
            f21 = sfe_dx[(j-1)*Nx+i];
            f22 = sfe_dx[j*Nx+i];
          }

          opdx = 1. / ((x2-x1)*(y2-y1)) * (
            f11 * (x2-opx) * (y2-opy) +
            f21 * (opx-x1) * (y2-opy) +
            f12 * (x2-opx) * (opy-y1) +
            f22 * (opx-x1) * (opy-y1));


          opd_pfdd[s] = opdx;
          index = ia*N2*N3 + ja*N3 + ka + s*lN1*N2*N3;
          dE_core[index] = opd_pfdd[s];

          }
        }
      }//ka
    }//s
}
    /* -----------------------------------------------------------------------
    core energy extended dislocations
    ---------------------------------------------------------------------*/
    /* Claire added extended energy 07/30/18 to account for partial dislocations
    and the Schoeck parameterization. Can be used for both fcc {111} planes and
    hcp {0001} basal planes to model partial dislocations*/
    void FFTW_Slab_Conduct::core_energy_extended()
    {
      int i, j, k, isa, index, index1, index2, index3, plane, num;
      int tag;
      //  int counter, marker, tag, countSF, cSF, count;

      int ND = dimension;
      int N1 = nx;
      double size = static_cast<double>(N1);
      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NP = num_planes;
      int NS = slip_systems;
      double dslip = material->dslip;

      int NT = solve->max_iter;
      int NSI = app->stoptime;

      double c0 = material->c0;
      double c1 = material->c1;
      double c2 = material->c2;
      double c3 = material->c3;
      double c4 = material->c4;
      double a1 = material->a1;
      double a3 = material->a3;
      double isf = material->isf;
      double usf = material->usf;
      double An = material->An;
      double Cn = material->Cn;
      double b = material->b;

      double mu = material->mu;
      double ll = material->ll;
      double young = material->young;
      double xnu = material->xnu;

      c0 = c0/(mu);
      c1 = c1/(mu);
      c2 = c2/(mu);
      c3 = c3/(mu);
      c4 = c4/(mu);
      a1 = a1/(mu);
      a3 = a3/(mu);
      isf = isf/(mu*dslip*b);
      An = An/(mu*dslip*b);
      Cn = (usf - (isf/2.0))/(mu*dslip*b);

      //  counter = 0;
      //  countSF = 0;
      //  marker = 0;

      /*as of now the planes need to be in the correct order initially for this to work
      also make sure to call this subroutine before the FFTW functions are called so that data_fftw = xi[na0] or double check that this is true.
      subroutine calculates core energy for each time step, need to called every time step and initialize every time step.*/

      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          data_core[index].re = data_fftw[index].re;
          data_core[index].im = data_fftw[index].im;
        }/*ijk*/
      }/*isa*/

      for(plane=0;plane<NP;plane++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          /*no derivatives taken use to calculate E_core*/

          fcore[index] = (c0 + c1*(cos(2.0*M_PI*(data_core[index1].re-data_core[index2].re)) + cos(2.0*M_PI*(data_core[index2].re-data_core[index3].re)) + cos(2.0*M_PI*(data_core[index3].re-data_core[index1].re))) + c2*(cos(2.0*M_PI*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + cos(2.0*M_PI*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + cos(2.0*M_PI*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + c3*(cos(4.0*M_PI*(data_core[index1].re-data_core[index2].re)) + cos(4.0*M_PI*(data_core[index2].re-data_core[index3].re)) + cos(4.0*M_PI*(data_core[index3].re-data_core[index1].re))) + c4*(cos(2.0*M_PI*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + cos(2.0*M_PI*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + cos(2*M_PI*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + cos(2.0*M_PI*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + cos(2.0*M_PI*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + cos(2.0*M_PI*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(sin(2.0*M_PI*(data_core[index1].re-data_core[index2].re)) + sin(2.0*M_PI*(data_core[index2].re-data_core[index3].re)) + sin(2.0*M_PI*(data_core[index3].re-data_core[index1].re))) + a3*(sin(4.0*M_PI*(data_core[index1].re-data_core[index2].re)) + sin(4.0*M_PI*(data_core[index2].re-data_core[index3].re)) + sin(4.0*M_PI*(data_core[index3].re-data_core[index1].re))))/(dslip*b);

          MPI_Reduce(&fcore[index], &f_core[index], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

          //            E_core += (f_core[index]/nsize);

          /*partial derivative wrt phase field 1*/
          df1core[index] = ((2.0*M_PI)/(dslip*b))*(c1*(sin(2.0*M_PI*(data_core[index2].re-data_core[index1].re)) + sin(2.0*M_PI*(data_core[index3].re-data_core[index1].re))) + c2*(2.0*sin(2.0*M_PI*(data_core[index2].re+data_core[index3].re-2.0*data_core[index1].re)) + sin(2.0*M_PI*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + sin(2.0*M_PI*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + (2.0*c3)*(sin(4.0*M_PI*(data_core[index2].re-data_core[index1].re)) + sin(4.0*M_PI*(data_core[index3].re-data_core[index1].re))) + c4*(3.0*sin(2.0*M_PI*(data_core[index2].re+2.0*data_core[index3].re-3.0*data_core[index1].re)) + 3.0*sin(2.0*M_PI*(2.0*data_core[index2].re+data_core[index3].re-3.0*data_core[index1].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + sin(2.0*M_PI*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + sin(2.0*M_PI*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(cos(2.0*M_PI*(data_core[index1].re-data_core[index2].re)) - cos(2.0*M_PI*(data_core[index3].re-data_core[index1].re))) + (2.0*a3)*(cos(4.0*M_PI*(data_core[index1].re-data_core[index2].re)) - cos(4.0*M_PI*(data_core[index3].re-data_core[index1].re))));

          /*partial derivative wrt phase field 2*/
          df2core[index] = ((2.0*M_PI)/(dslip*b))*(c1*(sin(2.0*M_PI*(data_core[index1].re-data_core[index2].re)) + sin(2.0*M_PI*(data_core[index3].re-data_core[index2].re))) + c2*(sin(2.0*M_PI*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + 2.0*sin(2.0*M_PI*(data_core[index3].re+data_core[index1].re-2.0*data_core[index2].re)) + sin(2.0*M_PI*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + (2.0*c3)*(sin(4.0*M_PI*(data_core[index1].re-data_core[index2].re)) + sin(4.0*M_PI*(data_core[index3].re-data_core[index2].re))) + c4*(sin(2.0*M_PI*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + 3.0*sin(2.0*M_PI*(data_core[index3].re+2.0*data_core[index1].re-3.0*data_core[index2].re)) + 3.0*sin(2.0*M_PI*(2.0*data_core[index3].re+data_core[index1].re-3.0*data_core[index2].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + sin(2.0*M_PI*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(cos(2.0*M_PI*(data_core[index2].re-data_core[index3].re)) - cos(2.0*M_PI*(data_core[index1].re-data_core[index2].re))) + (2.0*a3)*(cos(4.0*M_PI*(data_core[index2].re-data_core[index3].re)) - cos(4.0*M_PI*(data_core[index1].re-data_core[index2].re))));

          /*partial derivative wrt phase field 3*/
          df3core[index] = ((2.0*M_PI)/(dslip*b))*(c1*(sin(2.0*M_PI*(data_core[index2].re-data_core[index3].re)) + sin(2.0*M_PI*(data_core[index1].re-data_core[index3].re))) + c2*(sin(2.0*M_PI*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + sin(2.0*M_PI*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + 2.0*sin(2.0*M_PI*(data_core[index1].re+data_core[index2].re-2.0*data_core[index3].re))) + (2.0*c3)*(sin(4.0*M_PI*(data_core[index2].re-data_core[index3].re)) + sin(4.0*M_PI*(data_core[index1].re-data_core[index3].re))) + c4*(2.0*sin(2.0*M_PI*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + sin(2.0*M_PI*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + sin(2.0*M_PI*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + 3.0*sin(2.0*M_PI*(data_core[index1].re+2.0*data_core[index2].re-3.0*data_core[index3].re)) + 3.0*sin(2.0*M_PI*(2.0*data_core[index1].re+data_core[index2].re-3.0*data_core[index3].re))) + a1*(cos(2.0*M_PI*(data_core[index3].re-data_core[index1].re)) - cos(2.0*M_PI*(data_core[index2].re-data_core[index3].re))) + (2.0*a3)*(cos(4.0*M_PI*(data_core[index3].re-data_core[index1].re)) - cos(4.0*M_PI*(data_core[index2].re-data_core[index3].re))));
        }/*ijk*/
      }/*plane*/

      for(plane=0;plane<NP;plane++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          dE_core[index1] = df1core[index];
          dE_core[index2] = df2core[index];
          dE_core[index3] = df3core[index];
        }/*ijk*/
      }/*plane*/
      MPI_Barrier(MPI_COMM_WORLD);

    }

/* -----------------------------------------------------------------------
core energy extended dislocations (11 parameter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended11()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double c0 = material->c0;
  double c1 = material->c1;
  double c2 = material->c2;
  double c3 = material->c3;
  double c4 = material->c4;
  double c5 = material->c5;
  double c6 = material->c6;
  double c7 = material->c7;
  double c8 = material->c8;
  double c9 = material->c9;
  double c10 = material->c10;
  double c11 = 0.0;
  double b = material->b;

  double mu = material->mu;

  c0 = c0/(mu*dslip*b);
  c1 = c1/(mu*dslip*b);
  c2 = c2/(mu*dslip*b);
  c3 = c3/(mu*dslip*b);
  c4 = c4/(mu*dslip*b);
  c5 = c5/(mu*dslip*b);
  c6 = c6/(mu*dslip*b);
  c7 = c7/(mu*dslip*b);
  c8 = c8/(mu*dslip*b);
  c9 = c9/(mu*dslip*b);
  c10 = c10/(mu*dslip*b);
  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = c0 + c1 * (cos(2*M_PI*(xi2-xi3)) + cos(2*M_PI*(xi3-xi1)) + cos(2*M_PI*(xi1-xi2))) + c2 * (cos(2*M_PI*(2*xi1-xi2-xi3)) + cos(2*M_PI*(2*xi3-xi1-xi2)) + cos(2*M_PI*(2*xi2-xi1-xi3))) + c3 * (cos(4*M_PI*(xi2-xi3)) + cos(4*M_PI*(xi3-xi1)) + cos(4*M_PI*(xi1-xi2))) + c4 * (cos(2*M_PI*(3*xi1-xi2-2*xi3)) + cos(2*M_PI*(3*xi1-2*xi2-xi3)) + cos(2*M_PI*(2*xi1+xi2-3*xi3)) + cos(2*M_PI*(2*xi1-3*xi2+xi3)) + cos(2*M_PI*(xi1+2*xi2-3*xi3)) + cos(2*M_PI*(xi1-3*xi2+2*xi3))) + c5 * (sin(2*M_PI*(xi2-xi3)) + sin(2*M_PI*(xi3-xi1)) + sin(2*M_PI*(xi1-xi2))) + c6 * (sin(4*M_PI*(xi2-xi3)) + sin(4*M_PI*(xi3-xi1)) + sin(4*M_PI*(xi1-xi2))) + c7 * (cos(6*M_PI*(xi2-xi3)) + cos(6*M_PI*(xi3-xi1)) + cos(6*M_PI*(xi1-xi2))) + c8 * (cos(8*M_PI*(xi2-xi3)) + cos(8*M_PI*(xi3-xi1)) + cos(8*M_PI*(xi1-xi2))) + c9 * (cos(4*M_PI*(2*xi1-xi2-xi3)) + cos(4*M_PI*(xi1+xi2-2*xi3)) + cos(4*M_PI*(xi1-2*xi2+xi3))) + c10 *(cos(6*M_PI*(2*xi1-xi2-xi3)) + cos(6*M_PI*(xi1+xi2-2*xi3)) + cos(6*M_PI*(xi1-2*xi2+xi3)));

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = 2.0*M_PI*(c1 * (sin(2*M_PI*(xi2-xi1)) + sin(2*M_PI*(xi3-xi1))) + c2 * (-2*sin(2*M_PI*(2*xi1-xi2-xi3)) + sin(2*M_PI*(2*xi3-xi1-xi2)) + sin(2*M_PI*(2*xi2-xi1-xi3))) + c3 * (-2*sin(4*M_PI*(xi1-xi3)) - 2*sin(4*M_PI*(xi1-xi2))) + c4 * (-3*sin(2*M_PI*(3*xi1-xi2-2*xi3)) - 3*sin(2*M_PI*(3*xi1-2*xi2-xi3)) - 2*sin(2*M_PI*(2*xi1+xi2-3*xi3)) - 2*sin(2*M_PI*(2*xi1-3*xi2+xi3)) - sin(2*M_PI*(xi1+2*xi2-3*xi3)) + sin(2*M_PI*(-xi1+3*xi2-2*xi3))) + c5 * (-cos(2*M_PI*(xi3-xi1)) + cos(2*M_PI*(xi1-xi2))) + c6 * (-2*cos(4*M_PI*(xi3-xi1)) + 2*cos(4*M_PI*(xi1-xi2))) + c7 * (3*sin(6*M_PI*(xi3-xi1)) - 3*sin(6*M_PI*(xi1-xi2))) + c8 * (4*sin(8*M_PI*(xi3-xi1)) - 4*sin(8*M_PI*(xi1-xi2))) + c9 * (-4*sin(4*M_PI*(2*xi1-xi2-xi3)) - 2*sin(4*M_PI*(xi1+xi2-2*xi3)) + 2*sin(4*M_PI*(2*xi2-xi1-xi3))) + c10 *(-6*sin(6*M_PI*(2*xi1-xi2-xi3)) - 3*sin(6*M_PI*(xi1+xi2-2*xi3)) + 3*sin(6*M_PI*(2*xi2-xi1-xi3))));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = 2.0*M_PI*(c1 * (-sin(2*M_PI*(xi2-xi3)) - sin(2*M_PI*(xi2-xi1))) + c2 * (sin(2*M_PI*(2*xi1-xi2-xi3)) + sin(2*M_PI*(2*xi3-xi1-xi2)) - 2*sin(2*M_PI*(2*xi2-xi1-xi3))) + c3 * (-2*sin(4*M_PI*(xi2-xi3)) + 2*sin(4*M_PI*(xi1-xi2))) + c4 * (sin(2*M_PI*(3*xi1-xi2-2*xi3)) + 2*sin(2*M_PI*(3*xi1-2*xi2-xi3)) - sin(2*M_PI*(2*xi1+xi2-3*xi3)) + 3*sin(2*M_PI*(2*xi1-3*xi2+xi3)) - 2*sin(2*M_PI*(xi1+2*xi2-3*xi3)) - 3*sin(2*M_PI*(-xi1+3*xi2-2*xi3))) + c5 * (cos(2*M_PI*(xi2-xi3)) - cos(2*M_PI*(xi1-xi2))) + c6 * (2*cos(4*M_PI*(xi2-xi3)) - 2*cos(4*M_PI*(xi1-xi2))) + c7 * (-3*sin(6*M_PI*(xi2-xi3)) + 3*sin(6*M_PI*(xi1-xi2))) + c8 * (-4*sin(8*M_PI*(xi2-xi3)) + 4*sin(8*M_PI*(xi1-xi2))) + c9 * (2*sin(4*M_PI*(2*xi1-xi2-xi3)) - 2*sin(4*M_PI*(xi1+xi2-2*xi3)) - 4*sin(4*M_PI*(2*xi2-xi1-xi3))) + c10 *(3*sin(6*M_PI*(2*xi1-xi2-xi3)) - 3*sin(6*M_PI*(xi1+xi2-2*xi3)) - 6*sin(6*M_PI*(2*xi2-xi1-xi3))));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = 2.0*M_PI*(c1 * (sin(2*M_PI*(xi2-xi3)) + sin(2*M_PI*(xi1-xi3))) + c2 * (sin(2*M_PI*(2*xi1-xi2-xi3)) - 2*sin(2*M_PI*(2*xi3-xi1-xi2)) + sin(2*M_PI*(2*xi2-xi1-xi3))) + c3 * (2*sin(4*M_PI*(xi2-xi3)) + 2*sin(4*M_PI*(xi1-xi3))) + c4 * (2*sin(2*M_PI*(3*xi1-xi2-2*xi3)) + sin(2*M_PI*(3*xi1-2*xi2-xi3)) + 3*sin(2*M_PI*(2*xi1+xi2-3*xi3)) - sin(2*M_PI*(2*xi1-3*xi2+xi3)) + 3*sin(2*M_PI*(xi1+2*xi2-3*xi3)) + 2*sin(2*M_PI*(-xi1+3*xi2-2*xi3))) + c5 * (-cos(2*M_PI*(xi2-xi3)) + cos(2*M_PI*(xi3-xi1))) + c6 * (-2*cos(4*M_PI*(xi2-xi3)) + 2*cos(4*M_PI*(xi3-xi1))) + c7 * (3*sin(6*M_PI*(xi2-xi3)) - 3*sin(6*M_PI*(xi3-xi1))) + c8 * (4*sin(8*M_PI*(xi2-xi3)) - 4*sin(8*M_PI*(xi3-xi1))) + c9 * (2*sin(4*M_PI*(2*xi1-xi2-xi3)) + 4*sin(4*M_PI*(xi1+xi2-2*xi3)) + 2*sin(4*M_PI*(2*xi2-xi1-xi3))) + c10 *(3*sin(6*M_PI*(2*xi1-xi2-xi3)) + 6*sin(6*M_PI*(xi1+xi2-2*xi3)) + 3*sin(6*M_PI*(2*xi2-xi1-xi3))));

        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Ag (Multishifter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Ag()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  double omul = material->omul;


  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = -5.731787*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 93.85783*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 1.911371*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 97.20805*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 97.2037*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.909977*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.280206*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.175411*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.861596*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.933037*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.927822*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.934661*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.153347*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 57.00842*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 15.54031*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 1.002806*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 58.49648*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 58.65317*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.159734*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.404015*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.22774*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 16.0704*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 15.95516*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.236685*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.170065*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.210744*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.301471*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 220.1973;

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(62.16124*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 6.016836*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 116.99296*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 117.30634*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.958404*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 13.61606*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.91096*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 32.1408*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 31.91032*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.94674*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.68026*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.421488*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.602942*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 11.468226*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 194.4161*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 194.4074*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.459862*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.120824*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 24.701644*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.446384*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.732148*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.855644*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.869322*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(12.613388*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 114.01684*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 31.08062*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 2.005612*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 117.30634*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.638936*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 12.91096*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 32.1408*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 63.82064*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.47337*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.02039*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.842976*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.808826*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 22.927148*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 187.71566*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 3.822742*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 194.4074*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.639908*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 24.701644*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.723192*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 11.598222*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.711288*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.607966*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-12.613388*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 114.01684*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 31.08062*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 4.011224*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 116.99296*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.319468*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 13.61606*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 64.2816*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 31.91032*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.42011*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.34013*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.264464*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.205884*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 22.927148*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 187.71566*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 7.645484*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 194.4161*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.819954*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 25.120824*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.169576*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.866074*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 11.566932*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.738644*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          if(/*(lxs+i) == 0 && */((k < int(N3/8.0)) || k > int(N3*7.0/8.0))){
            dE_core[index1] = omul*dE_core[index1];
            dE_core[index2] = omul*dE_core[index2];
            dE_core[index3] = omul*dE_core[index3];
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Au (Multishifter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Au()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  double omul = material->omul;

  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = -7.671776*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 73.43161*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 1.921959*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 71.45345*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 71.42202*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.935711*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.606874*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.653461*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.052228*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.053703*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.081598*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.067353*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.584475*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 46.37485*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 12.85583*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 1.346036*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 45.58045*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 45.60515*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.414934*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.621059*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.598583*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 13.71983*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 13.68223*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.400701*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.451948*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.420884*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.390111*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 176.3105;

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(51.42332*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 8.076216*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 91.1609*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 91.2103*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.489604*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.484236*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.394332*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 27.43966*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 27.36446*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.602804*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.807792*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.841768*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.780222*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 11.531754*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 142.9069*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 142.84404*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.614266*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 30.427496*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 30.613844*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.208912*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.214812*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.163196*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.134706*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(14.3379*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 92.7497*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 25.71166*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 2.692072*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 91.2103*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.659736*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.394332*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 27.43966*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 54.72892*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.801402*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.711688*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.683536*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.340666*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 30.687104*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 146.86322*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 3.843918*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 142.84404*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.742844*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 30.613844*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.104456*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.322218*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.326392*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 12.404118*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-14.3379*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 92.7497*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 25.71166*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 5.384144*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 91.1609*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.829868*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.484236*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 54.87932*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 27.36446*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.404206*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.903896*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.525304*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.560444*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 30.687104*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 146.86322*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 7.687836*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 142.9069*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.871422*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 30.427496*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 12.313368*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.107406*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.489588*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.269412*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          if(/*(lxs+i) == 0 && */((k < int(N3/8.0)) || k > int(N3*7.0/8.0))){
            dE_core[index1] = omul*dE_core[index1];
            dE_core[index2] = omul*dE_core[index2];
            dE_core[index3] = omul*dE_core[index3];
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Cu (Multishifter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Cu()
{
  int i, j, k, isa, index, index1, index2, index3, ptindex, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  double omul = material->omul;

  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;
          ptindex = i*N2*N3 + j*N3 + k;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = -8.329971*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 165.3755*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 2.411782*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 170.1183*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 170.1881*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.175776*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.01989*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.041115*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.13134*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.295692*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.344248*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.430872*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.074517*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 104.3822*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 22.1248*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 1.33398*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 105.7677*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 105.8889*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.293769*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.398633*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.377277*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 22.91647*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 22.74947*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.467355*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.50688*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.278417*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.406927*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 379.3413;
          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(88.4992*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 8.00388*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 211.5354*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 211.7778*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.762614*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 17.594532*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 17.509108*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 45.83294*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 45.49894*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.86942*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.02752*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.556834*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.813854*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.470692*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 340.2366*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 340.3762*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 13.054656*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 36.07956*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 36.16446*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.52536*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.182768*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.688496*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.861744*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(16.298068*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 208.7644*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 44.2496*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 2.66796*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 211.7778*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.175076*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 17.509108*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 45.83294*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 90.99788*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.93471*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.04128*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.113668*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.441562*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 33.319884*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 330.751*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 4.823564*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 340.3762*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.703104*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 36.16446*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.26268*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 13.774152*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.376992*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.585232*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-16.298068*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 208.7644*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 44.2496*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 5.33592*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 211.5354*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.587538*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 17.594532*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 91.66588*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 45.49894*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.80413*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.01376*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.670502*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.627708*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 33.319884*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 330.751*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 9.647128*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 340.2366*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.351552*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 36.07956*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 12.78804*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.591384*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.065488*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.723488*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          
          if(/*(lxs+i) == 0 && */((k < int(N3/8.0)) || k > int(N3*7.0/8.0))){
          // if(/*(lxs+i) <= int(N1/4.0) && */((k < int(N3*3.0/8.0) && k >= int(N3*1.0/4.0)) || (k > int(N3*5.0/8.0) && k <= int(N3*3.0/4.0)))){ // 64x64x64 --> 64x64x128
            fcore[index] = omul*fcore[index];
            dE_core[index1] = omul*dE_core[index1];
            dE_core[index2] = omul*dE_core[index2];
            dE_core[index3] = omul*dE_core[index3];
          }
          fcore[index] = fcore[index]*0.001;
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Cu (Multishifter) with Twin Boundary testing (GSFE & GTFE surface calculated by Ian Chesser from Mishin Cu EAM1, Multishifter series)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Cu_gbtest()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  
  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = 1.014712*sin(M_PI*(-10.0*xi2 + 10.0*xi3)) + 1.141076*sin(M_PI*(-9.0*xi2 + 9.0*xi3)) + 1.067101*sin(M_PI*(-8.0*xi2 + 8.0*xi3)) + 1.516156*sin(M_PI*(-7.0*xi2 + 7.0*xi3)) + 1.999785*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) + 2.525744*sin(M_PI*(-5.0*xi2 + 5.0*xi3)) - 7.503378*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 4.994296*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) - 131.3614*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.38802*sin(M_PI*(-1.0*xi2 + xi3)) + 1.062852*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 7.493571*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 3.418324*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 135.3559*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 161.6725*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.182386*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.251403*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.609219*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.92338*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.102621*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.062878*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.740584*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.890952*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.746084*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.736207*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.210922*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.440973*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.413646*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.04073*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.434426*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.703986*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) + 3.906802*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 1.336228*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) - 106.0941*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 12.31535*cos(M_PI*(-1.0*xi2 + xi3)) - 7.048628*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 24.28289*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 1.969479*cos(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 1.447656*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.185268*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 109.657*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 64.86831*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.531094*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.277132*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.036736*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.113594*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.61421*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.0329*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.046529*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.446845*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.574019*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 343.841;

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(14.097256*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 97.13156*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 11.816874*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 11.581248*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.741072*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 219.314*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 129.73662*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.124376*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.108528*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 32.146944*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.681564*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.22842*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 50.0658*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.186116*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.89369*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.148038*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.125704*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 29.974284*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 13.673296*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 270.7118*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 323.345*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 20.729544*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.508418*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.873752*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 43.69352*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.205242*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.125756*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 30.962336*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.781904*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.492168*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.944828*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.421844*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.881946*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 13.654584*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.08146*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.868852*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(10.223916*M_PI*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) + 15.627208*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 4.008684*M_PI*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) - 212.1882*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 12.31535*M_PI*sin(M_PI*(-1.0*xi2 + xi3)) - 7.048628*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 48.56578*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 5.908437*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) - 4.342968*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.185268*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 129.73662*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.593282*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 32.146944*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.22842*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 100.1316*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.093058*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.78738*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 15.444114*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.14712*M_PI*cos(M_PI*(-10.0*xi2 + 10.0*xi3)) - 10.269684*M_PI*cos(M_PI*(-9.0*xi2 + 9.0*xi3)) - 8.536808*M_PI*cos(M_PI*(-8.0*xi2 + 8.0*xi3)) - 10.613092*M_PI*cos(M_PI*(-7.0*xi2 + 7.0*xi3)) - 11.99871*M_PI*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) - 12.62872*M_PI*cos(M_PI*(-5.0*xi2 + 5.0*xi3)) + 30.013512*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 14.982888*M_PI*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) + 262.7228*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 11.38802*M_PI*cos(M_PI*(-1.0*xi2 + xi3)) - 1.062852*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 14.987142*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 3.418324*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 323.345*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.547158*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.005612*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.046095*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.102621*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.188634*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 30.962336*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.781904*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 22.984336*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.472414*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.632766*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.204865*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 20.481876*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.16292*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 20.606556*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-10.223916*M_PI*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) - 15.627208*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 4.008684*M_PI*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) + 212.1882*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 12.31535*M_PI*sin(M_PI*(-1.0*xi2 + xi3)) - 7.048628*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 48.56578*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 5.908437*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) - 7.23828*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.555804*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 219.314*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.531094*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.108528*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.681564*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 102.45684*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 50.0658*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 12.279174*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.68107*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.296076*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.14712*M_PI*cos(M_PI*(-10.0*xi2 + 10.0*xi3)) + 10.269684*M_PI*cos(M_PI*(-9.0*xi2 + 9.0*xi3)) + 8.536808*M_PI*cos(M_PI*(-8.0*xi2 + 8.0*xi3)) + 10.613092*M_PI*cos(M_PI*(-7.0*xi2 + 7.0*xi3)) + 11.99871*M_PI*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) + 12.62872*M_PI*cos(M_PI*(-5.0*xi2 + 5.0*xi3)) - 30.013512*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 14.982888*M_PI*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) - 262.7228*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.38802*M_PI*cos(M_PI*(-1.0*xi2 + xi3)) - 1.062852*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 14.987142*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 10.254972*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 270.7118*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.182386*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.502806*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.827657*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 43.69352*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.307863*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.062878*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.563808*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.492168*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.417242*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.05461*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.322919*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.827292*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 24.24438*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 13.737704*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          // if(/*(lxs+i) == 0 &&*/ ((k == int(N3/8.0)) || k == int(N3*7.0/8.0))){
          // fcore[index] = 0.2815984*sin(M_PI*(-21.0*xi2 + 21.0*xi3)) + 0.3195724*sin(M_PI*(-20.0*xi2 + 20.0*xi3)) + 0.3576857*sin(M_PI*(-19.0*xi2 + 19.0*xi3)) + 0.3986237*sin(M_PI*(-18.0*xi2 + 18.0*xi3)) + 0.4450561*sin(M_PI*(-17.0*xi2 + 17.0*xi3)) + 0.4915005*sin(M_PI*(-16.0*xi2 + 16.0*xi3)) + 0.5500461*sin(M_PI*(-15.0*xi2 + 15.0*xi3)) + 0.6084532*sin(M_PI*(-14.0*xi2 + 14.0*xi3)) + 0.6845937*sin(M_PI*(-13.0*xi2 + 13.0*xi3)) + 0.7486672*sin(M_PI*(-12.0*xi2 + 12.0*xi3)) + 0.8539095*sin(M_PI*(-11.0*xi2 + 11.0*xi3)) + 1.043075*sin(M_PI*(-10.0*xi2 + 10.0*xi3)) + 1.16352*sin(M_PI*(-9.0*xi2 + 9.0*xi3)) + 1.175222*sin(M_PI*(-8.0*xi2 + 8.0*xi3)) + 1.581174*sin(M_PI*(-7.0*xi2 + 7.0*xi3)) + 1.780717*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) + 2.563205*sin(M_PI*(-5.0*xi2 + 5.0*xi3)) - 7.011371*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 5.230236*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) - 136.0463*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.62517*sin(M_PI*(-1.0*xi2 + xi3)) + 1.029216*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 7.493423*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 0.9430937*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) - 0.3373575*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 0.3368915*sin(M_PI*(10*xi1 - 5.0*xi2 - 5.0*xi3)) - 0.2705567*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.7642227*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.7785793*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.634125*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 140.0523*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 164.1832*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.255622*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.284622*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.601067*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.7705179*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.4959997*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 12.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.2802797*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 14.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.9057483*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.51417*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.261301*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.226851*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.709637*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.4928758*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.3274968*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.901482*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.754829*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.4688069*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.6060518*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.701923*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.252809*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.464715*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.464299*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.012999*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.363912*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.2583958*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.6902584*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.7786087*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.8244815*sin(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.8581661*sin(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.5887556*sin(16.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.5828665*sin(16.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.5430441*sin(18.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.5268043*sin(18.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.4855487*sin(20.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.4494857*sin(20.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.4321445*sin(22.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.3899314*sin(22.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.3987247*sin(24.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.3495924*sin(24.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.3614754*sin(26.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.3072719*sin(26.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.3310763*sin(28.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.2738114*sin(28.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.3029041*sin(30.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.2776648*sin(32.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.2552034*sin(34.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.5489964*cos(M_PI*(-30.0*xi2 + 30.0*xi3)) + 0.5487554*cos(M_PI*(-29.0*xi2 + 29.0*xi3)) + 0.549885*cos(M_PI*(-28.0*xi2 + 28.0*xi3)) + 0.550238*cos(M_PI*(-27.0*xi2 + 27.0*xi3)) + 0.5509822*cos(M_PI*(-26.0*xi2 + 26.0*xi3)) + 0.55131*cos(M_PI*(-25.0*xi2 + 25.0*xi3)) + 0.5537728*cos(M_PI*(-24.0*xi2 + 24.0*xi3)) + 0.5538111*cos(M_PI*(-23.0*xi2 + 23.0*xi3)) + 0.5588332*cos(M_PI*(-22.0*xi2 + 22.0*xi3)) + 0.5586842*cos(M_PI*(-21.0*xi2 + 21.0*xi3)) + 0.5641599*cos(M_PI*(-20.0*xi2 + 20.0*xi3)) + 0.5655853*cos(M_PI*(-19.0*xi2 + 19.0*xi3)) + 0.5694084*cos(M_PI*(-18.0*xi2 + 18.0*xi3)) + 0.5739412*cos(M_PI*(-17.0*xi2 + 17.0*xi3)) + 0.5791474*cos(M_PI*(-16.0*xi2 + 16.0*xi3)) + 0.5859943*cos(M_PI*(-15.0*xi2 + 15.0*xi3)) + 0.5855692*cos(M_PI*(-14.0*xi2 + 14.0*xi3)) + 0.5980104*cos(M_PI*(-13.0*xi2 + 13.0*xi3)) + 0.6393258*cos(M_PI*(-12.0*xi2 + 12.0*xi3)) + 0.6391327*cos(M_PI*(-11.0*xi2 + 11.0*xi3)) + 0.6662425*cos(M_PI*(-10.0*xi2 + 10.0*xi3)) + 0.6976632*cos(M_PI*(-9.0*xi2 + 9.0*xi3)) + 0.4962071*cos(M_PI*(-8.0*xi2 + 8.0*xi3)) + 0.7116295*cos(M_PI*(-7.0*xi2 + 7.0*xi3)) + 1.569524*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) + 0.8473332*cos(M_PI*(-5.0*xi2 + 5.0*xi3)) + 5.060236*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 1.248211*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) - 99.42641*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 11.31079*cos(M_PI*(-1.0*xi2 + xi3)) - 6.82559*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 24.28249*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 1.957386*cos(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 0.4944981*cos(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 0.3628136*cos(M_PI*(10*xi1 - 5.0*xi2 - 5.0*xi3)) + 0.3324455*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 14.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.5326072*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 12.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.7525676*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.420141*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.022186*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.998882*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 102.9532*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 57.07792*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.255523*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.50673*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.4526047*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.821895*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.7653937*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.035899*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.7444321*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.67805*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.09701*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.5813453*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.9016895*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.186672*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.907543*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.4977697*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.2838154*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.596127*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.708482*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.2663616*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.4037039*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.3522789*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 344.0497;

          // /*partial derivative wrt phase field 1*/
          // dE_core[index1] = (0.001/(mu*dslip*b))*(13.65118*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 97.12996*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 11.744316*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) - 3.9559848*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 3.628136*M_PI*sin(M_PI*(10*xi1 - 5.0*xi2 - 5.0*xi3)) + 4.654237*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 14.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.3912864*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 12.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.525676*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.361128*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.133116*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.995528*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 205.9064*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 114.15584*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.022092*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.02692*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.9052094*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 35.28758*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.5923622*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.215394*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.9777284*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.3561*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 50.19402*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.3253812*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.410137*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.746688*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.815086*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.9955394*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.1352616*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.192254*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.416964*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.0654464*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.8074078*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.4091156*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.058432*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 29.973692*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 5.6585622*M_PI*cos(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) - 2.69886*M_PI*cos(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 3.368915*M_PI*cos(M_PI*(10*xi1 - 5.0*xi2 - 5.0*xi3)) + 2.705567*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.1137816*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.6714758*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.5365*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 280.1046*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 328.3664*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 21.022488*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.707732*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.808536*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.705179*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.9519964*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 12.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.9239158*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 14.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.4344898*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 42.05668*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.522602*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.453702*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 26.838548*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.9572548*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.9649808*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.802964*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.509658*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.8752276*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.6363108*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.807692*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.505618*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.92943*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 13.857196*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.025998*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.727824*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.0335832*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.3805168*M_PI*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.5572174*M_PI*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.648963*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.7163322*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.1775112*M_PI*cos(16.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.165733*M_PI*cos(16.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.0860882*M_PI*cos(18.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.0536086*M_PI*cos(18.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.9710974*M_PI*cos(20.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.8989714*M_PI*cos(20.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.864289*M_PI*cos(22.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.7798628*M_PI*cos(22.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.7974494*M_PI*cos(24.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.6991848*M_PI*cos(24.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.7229508*M_PI*cos(26.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.6145438*M_PI*cos(26.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.6621526*M_PI*cos(28.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.5476228*M_PI*cos(28.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.6058082*M_PI*cos(30.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.5553296*M_PI*cos(32.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.5104068*M_PI*cos(34.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          // /*partial derivative wrt phase field 2*/
          // dE_core[index2] = (0.001/(mu*dslip*b))*(16.469892*M_PI*sin(M_PI*(-30.0*xi2 + 30.0*xi3)) + 15.9139066*M_PI*sin(M_PI*(-29.0*xi2 + 29.0*xi3)) + 15.39678*M_PI*sin(M_PI*(-28.0*xi2 + 28.0*xi3)) + 14.856426*M_PI*sin(M_PI*(-27.0*xi2 + 27.0*xi3)) + 14.3255372*M_PI*sin(M_PI*(-26.0*xi2 + 26.0*xi3)) + 13.78275*M_PI*sin(M_PI*(-25.0*xi2 + 25.0*xi3)) + 13.2905472*M_PI*sin(M_PI*(-24.0*xi2 + 24.0*xi3)) + 12.7376553*M_PI*sin(M_PI*(-23.0*xi2 + 23.0*xi3)) + 12.2943304*M_PI*sin(M_PI*(-22.0*xi2 + 22.0*xi3)) + 11.7323682*M_PI*sin(M_PI*(-21.0*xi2 + 21.0*xi3)) + 11.283198*M_PI*sin(M_PI*(-20.0*xi2 + 20.0*xi3)) + 10.7461207*M_PI*sin(M_PI*(-19.0*xi2 + 19.0*xi3)) + 10.2493512*M_PI*sin(M_PI*(-18.0*xi2 + 18.0*xi3)) + 9.7570004*M_PI*sin(M_PI*(-17.0*xi2 + 17.0*xi3)) + 9.2663584*M_PI*sin(M_PI*(-16.0*xi2 + 16.0*xi3)) + 8.7899145*M_PI*sin(M_PI*(-15.0*xi2 + 15.0*xi3)) + 8.1979688*M_PI*sin(M_PI*(-14.0*xi2 + 14.0*xi3)) + 7.7741352*M_PI*sin(M_PI*(-13.0*xi2 + 13.0*xi3)) + 7.6719096*M_PI*sin(M_PI*(-12.0*xi2 + 12.0*xi3)) + 7.0304597*M_PI*sin(M_PI*(-11.0*xi2 + 11.0*xi3)) + 6.662425*M_PI*sin(M_PI*(-10.0*xi2 + 10.0*xi3)) + 6.2789688*M_PI*sin(M_PI*(-9.0*xi2 + 9.0*xi3)) + 3.9696568*M_PI*sin(M_PI*(-8.0*xi2 + 8.0*xi3)) + 4.9814065*M_PI*sin(M_PI*(-7.0*xi2 + 7.0*xi3)) + 9.417144*M_PI*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) + 4.236666*M_PI*sin(M_PI*(-5.0*xi2 + 5.0*xi3)) + 20.240944*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 3.744633*M_PI*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) - 198.85282*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 11.31079*M_PI*sin(M_PI*(-1.0*xi2 + xi3)) - 6.82559*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 48.56498*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 5.872158*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 1.9779924*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 1.814068*M_PI*sin(M_PI*(10*xi1 - 5.0*xi2 - 5.0*xi3)) - 1.994673*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 14.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.663036*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 12.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.0102704*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.260423*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.044372*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.998882*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 114.15584*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.766569*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.3578141*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 35.28758*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.8269685*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.7444321*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.3561*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 100.38804*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.9067265*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.410137*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.373344*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.722629*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.4888485*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.7028924*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.384508*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 16.250892*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.0654464*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.0185195*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.8182312*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.9135664*M_PI*cos(M_PI*(-21.0*xi2 + 21.0*xi3)) - 6.391448*M_PI*cos(M_PI*(-20.0*xi2 + 20.0*xi3)) - 6.7960283*M_PI*cos(M_PI*(-19.0*xi2 + 19.0*xi3)) - 7.1752266*M_PI*cos(M_PI*(-18.0*xi2 + 18.0*xi3)) - 7.5659537*M_PI*cos(M_PI*(-17.0*xi2 + 17.0*xi3)) - 7.864008*M_PI*cos(M_PI*(-16.0*xi2 + 16.0*xi3)) - 8.2506915*M_PI*cos(M_PI*(-15.0*xi2 + 15.0*xi3)) - 8.5183448*M_PI*cos(M_PI*(-14.0*xi2 + 14.0*xi3)) - 8.8997181*M_PI*cos(M_PI*(-13.0*xi2 + 13.0*xi3)) - 8.9840064*M_PI*cos(M_PI*(-12.0*xi2 + 12.0*xi3)) - 9.3930045*M_PI*cos(M_PI*(-11.0*xi2 + 11.0*xi3)) - 10.43075*M_PI*cos(M_PI*(-10.0*xi2 + 10.0*xi3)) - 10.47168*M_PI*cos(M_PI*(-9.0*xi2 + 9.0*xi3)) - 9.401776*M_PI*cos(M_PI*(-8.0*xi2 + 8.0*xi3)) - 11.068218*M_PI*cos(M_PI*(-7.0*xi2 + 7.0*xi3)) - 10.684302*M_PI*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) - 12.816025*M_PI*cos(M_PI*(-5.0*xi2 + 5.0*xi3)) + 28.045484*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 15.690708*M_PI*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) + 272.0926*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 11.62517*M_PI*cos(M_PI*(-1.0*xi2 + xi3)) - 1.029216*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 14.986846*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 2.8292811*M_PI*cos(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 1.34943*M_PI*cos(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 1.6844575*M_PI*cos(M_PI*(10*xi1 - 5.0*xi2 - 5.0*xi3)) - 1.0822268*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.2926681*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.5571586*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.634125*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 328.3664*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.766866*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.138488*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.005335*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.6231074*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.4719979*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 12.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.2422376*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 14.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.9057483*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.261301*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.680553*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 26.838548*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.464379*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.802964*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 23.019316*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.3440345*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.6363108*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.403846*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.758427*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.323575*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 20.785794*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.051996*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 20.183472*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.8087706*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.451292*M_PI*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.4502609*M_PI*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.946889*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.8653288*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.1212892*M_PI*cos(16.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.2457985*M_PI*cos(16.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.3443528*M_PI*cos(18.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.268043*M_PI*cos(18.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.3699383*M_PI*cos(20.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.9443427*M_PI*cos(20.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.321445*M_PI*cos(22.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.6791768*M_PI*cos(22.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.3859717*M_PI*cos(24.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.5447012*M_PI*cos(24.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.3377048*M_PI*cos(26.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.3018066*M_PI*cos(26.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.3039919*M_PI*cos(28.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.107171*M_PI*cos(28.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.2406574*M_PI*cos(30.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.164972*M_PI*cos(32.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.0832544*M_PI*cos(34.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          // /*partial derivative wrt phase field 3*/
          // dE_core[index3] = (0.001/(mu*dslip*b))*(-16.469892*M_PI*sin(M_PI*(-30.0*xi2 + 30.0*xi3)) - 15.9139066*M_PI*sin(M_PI*(-29.0*xi2 + 29.0*xi3)) - 15.39678*M_PI*sin(M_PI*(-28.0*xi2 + 28.0*xi3)) - 14.856426*M_PI*sin(M_PI*(-27.0*xi2 + 27.0*xi3)) - 14.3255372*M_PI*sin(M_PI*(-26.0*xi2 + 26.0*xi3)) - 13.78275*M_PI*sin(M_PI*(-25.0*xi2 + 25.0*xi3)) - 13.2905472*M_PI*sin(M_PI*(-24.0*xi2 + 24.0*xi3)) - 12.7376553*M_PI*sin(M_PI*(-23.0*xi2 + 23.0*xi3)) - 12.2943304*M_PI*sin(M_PI*(-22.0*xi2 + 22.0*xi3)) - 11.7323682*M_PI*sin(M_PI*(-21.0*xi2 + 21.0*xi3)) - 11.283198*M_PI*sin(M_PI*(-20.0*xi2 + 20.0*xi3)) - 10.7461207*M_PI*sin(M_PI*(-19.0*xi2 + 19.0*xi3)) - 10.2493512*M_PI*sin(M_PI*(-18.0*xi2 + 18.0*xi3)) - 9.7570004*M_PI*sin(M_PI*(-17.0*xi2 + 17.0*xi3)) - 9.2663584*M_PI*sin(M_PI*(-16.0*xi2 + 16.0*xi3)) - 8.7899145*M_PI*sin(M_PI*(-15.0*xi2 + 15.0*xi3)) - 8.1979688*M_PI*sin(M_PI*(-14.0*xi2 + 14.0*xi3)) - 7.7741352*M_PI*sin(M_PI*(-13.0*xi2 + 13.0*xi3)) - 7.6719096*M_PI*sin(M_PI*(-12.0*xi2 + 12.0*xi3)) - 7.0304597*M_PI*sin(M_PI*(-11.0*xi2 + 11.0*xi3)) - 6.662425*M_PI*sin(M_PI*(-10.0*xi2 + 10.0*xi3)) - 6.2789688*M_PI*sin(M_PI*(-9.0*xi2 + 9.0*xi3)) - 3.9696568*M_PI*sin(M_PI*(-8.0*xi2 + 8.0*xi3)) - 4.9814065*M_PI*sin(M_PI*(-7.0*xi2 + 7.0*xi3)) - 9.417144*M_PI*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) - 4.236666*M_PI*sin(M_PI*(-5.0*xi2 + 5.0*xi3)) - 20.240944*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 3.744633*M_PI*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) + 198.85282*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.31079*M_PI*sin(M_PI*(-1.0*xi2 + xi3)) - 6.82559*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 48.56498*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 5.872158*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 1.9779924*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 1.814068*M_PI*sin(M_PI*(10*xi1 - 5.0*xi2 - 5.0*xi3)) - 2.659564*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 14.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.7282504*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 12.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.5154056*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.100705*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.088744*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 11.996646*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 205.9064*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.255523*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.02692*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.4526047*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.7653937*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.215394*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.7221605*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 102.7122*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 50.19402*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.5813453*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 13.120032*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.537715*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.4933091*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.5676308*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.576762*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.833928*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.1308928*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.8259273*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.4091156*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.9135664*M_PI*cos(M_PI*(-21.0*xi2 + 21.0*xi3)) + 6.391448*M_PI*cos(M_PI*(-20.0*xi2 + 20.0*xi3)) + 6.7960283*M_PI*cos(M_PI*(-19.0*xi2 + 19.0*xi3)) + 7.1752266*M_PI*cos(M_PI*(-18.0*xi2 + 18.0*xi3)) + 7.5659537*M_PI*cos(M_PI*(-17.0*xi2 + 17.0*xi3)) + 7.864008*M_PI*cos(M_PI*(-16.0*xi2 + 16.0*xi3)) + 8.2506915*M_PI*cos(M_PI*(-15.0*xi2 + 15.0*xi3)) + 8.5183448*M_PI*cos(M_PI*(-14.0*xi2 + 14.0*xi3)) + 8.8997181*M_PI*cos(M_PI*(-13.0*xi2 + 13.0*xi3)) + 8.9840064*M_PI*cos(M_PI*(-12.0*xi2 + 12.0*xi3)) + 9.3930045*M_PI*cos(M_PI*(-11.0*xi2 + 11.0*xi3)) + 10.43075*M_PI*cos(M_PI*(-10.0*xi2 + 10.0*xi3)) + 10.47168*M_PI*cos(M_PI*(-9.0*xi2 + 9.0*xi3)) + 9.401776*M_PI*cos(M_PI*(-8.0*xi2 + 8.0*xi3)) + 11.068218*M_PI*cos(M_PI*(-7.0*xi2 + 7.0*xi3)) + 10.684302*M_PI*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) + 12.816025*M_PI*cos(M_PI*(-5.0*xi2 + 5.0*xi3)) - 28.045484*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 15.690708*M_PI*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) - 272.0926*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.62517*M_PI*cos(M_PI*(-1.0*xi2 + xi3)) - 1.029216*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 14.986846*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 2.8292811*M_PI*cos(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 1.34943*M_PI*cos(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 1.6844575*M_PI*cos(M_PI*(10*xi1 - 5.0*xi2 - 5.0*xi3)) - 1.6233402*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.8211135*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.1143172*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.902375*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 280.1046*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.255622*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.569244*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.803201*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.0820716*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 10.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.4799985*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 12.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.6816782*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 14.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.5287415*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 42.05668*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.783903*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.226851*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 0.4928758*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.9649808*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.605928*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.509658*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.4688069*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.211538*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.264045*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.394145*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.928598*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 24.077994*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 13.455648*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 0.7751874*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.8318088*M_PI*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.8930435*M_PI*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.595852*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.1489966*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.2988004*M_PI*cos(16.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.0800655*M_PI*cos(16.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.430441*M_PI*cos(18.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.2144344*M_PI*cos(18.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.3410357*M_PI*cos(20.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.0453713*M_PI*cos(20.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.185734*M_PI*cos(22.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.899314*M_PI*cos(22.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.1834211*M_PI*cos(24.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.8455164*M_PI*cos(24.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.0606556*M_PI*cos(26.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.6872628*M_PI*cos(26.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.9661445*M_PI*cos(28.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.5595482*M_PI*cos(28.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.8464656*M_PI*cos(30.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.7203016*M_PI*cos(32.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.5936612*M_PI*cos(34.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          // }
          if(/*(lxs+i) == 0 && */((k == int(N3/8.0)) || k == int(N3*7.0/8.0))){ // LESS TERMS
          fcore[index] = 1.043075*sin(M_PI*(-10.0*xi2 + 10.0*xi3)) + 1.16352*sin(M_PI*(-9.0*xi2 + 9.0*xi3)) + 1.175222*sin(M_PI*(-8.0*xi2 + 8.0*xi3)) + 1.581174*sin(M_PI*(-7.0*xi2 + 7.0*xi3)) + 1.780717*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) + 2.563205*sin(M_PI*(-5.0*xi2 + 5.0*xi3)) - 7.011371*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 5.230236*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) - 136.0463*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.62517*sin(M_PI*(-1.0*xi2 + xi3)) + 1.029216*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 7.493423*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 3.634125*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 140.0523*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 164.1832*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.255622*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.284622*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.601067*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.51417*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.261301*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.226851*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.709637*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.901482*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.754829*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.701923*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.252809*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.464715*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.464299*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.012999*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.363912*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.569524*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) + 5.060236*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 1.248211*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) - 99.42641*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 11.31079*cos(M_PI*(-1.0*xi2 + xi3)) - 6.82559*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 24.28249*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 1.957386*cos(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 1.420141*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.022186*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.998882*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 102.9532*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 57.07792*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.255523*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.50673*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.821895*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.035899*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.67805*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.09701*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.186672*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.596127*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.708482*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 344.0497;

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(13.65118*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 97.12996*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 11.744316*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) + 11.361128*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.133116*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.995528*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 205.9064*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 114.15584*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.022092*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.02692*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 35.28758*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.215394*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.3561*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 50.19402*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.746688*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.192254*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.416964*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.058432*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) + 29.973692*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 14.5365*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 280.1046*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 328.3664*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 21.022488*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.707732*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.808536*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 42.05668*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.522602*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.453702*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 26.838548*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.802964*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.509658*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.807692*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.505618*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.92943*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 13.857196*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.025998*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.727824*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(9.417144*M_PI*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) + 20.240944*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 3.744633*M_PI*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) - 198.85282*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 11.31079*M_PI*sin(M_PI*(-1.0*xi2 + xi3)) - 6.82559*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 48.56498*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 5.872158*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) - 4.260423*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.044372*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.998882*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 114.15584*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.766569*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 35.28758*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.3561*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 100.38804*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.373344*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.384508*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 16.250892*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.43075*M_PI*cos(M_PI*(-10.0*xi2 + 10.0*xi3)) - 10.47168*M_PI*cos(M_PI*(-9.0*xi2 + 9.0*xi3)) - 9.401776*M_PI*cos(M_PI*(-8.0*xi2 + 8.0*xi3)) - 11.068218*M_PI*cos(M_PI*(-7.0*xi2 + 7.0*xi3)) - 10.684302*M_PI*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) - 12.816025*M_PI*cos(M_PI*(-5.0*xi2 + 5.0*xi3)) + 28.045484*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 15.690708*M_PI*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) + 272.0926*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 11.62517*M_PI*cos(M_PI*(-1.0*xi2 + xi3)) - 1.029216*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 14.986846*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 3.634125*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 328.3664*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.766866*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.138488*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.005335*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.261301*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.680553*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 26.838548*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.802964*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 23.019316*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.403846*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.758427*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.323575*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 20.785794*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.051996*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 20.183472*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-9.417144*M_PI*sin(M_PI*(-6.0*xi2 + 6.0*xi3)) - 20.240944*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 3.744633*M_PI*sin(M_PI*(-3.0*xi2 + 3.0*xi3)) + 198.85282*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.31079*M_PI*sin(M_PI*(-1.0*xi2 + xi3)) - 6.82559*M_PI*sin(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 48.56498*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 5.872158*M_PI*sin(M_PI*(6*xi1 - 3.0*xi2 - 3.0*xi3)) - 7.100705*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.088744*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 11.996646*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 205.9064*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.255523*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.02692*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.215394*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 102.7122*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 50.19402*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 13.120032*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.576762*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.833928*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.43075*M_PI*cos(M_PI*(-10.0*xi2 + 10.0*xi3)) + 10.47168*M_PI*cos(M_PI*(-9.0*xi2 + 9.0*xi3)) + 9.401776*M_PI*cos(M_PI*(-8.0*xi2 + 8.0*xi3)) + 11.068218*M_PI*cos(M_PI*(-7.0*xi2 + 7.0*xi3)) + 10.684302*M_PI*cos(M_PI*(-6.0*xi2 + 6.0*xi3)) + 12.816025*M_PI*cos(M_PI*(-5.0*xi2 + 5.0*xi3)) - 28.045484*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 15.690708*M_PI*cos(M_PI*(-3.0*xi2 + 3.0*xi3)) - 272.0926*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.62517*M_PI*cos(M_PI*(-1.0*xi2 + xi3)) - 1.029216*M_PI*cos(M_PI*(2*xi1 - 1.0*xi2 - 1.0*xi3)) - 14.986846*M_PI*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 10.902375*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 280.1046*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.255622*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.569244*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.803201*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 42.05668*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.783903*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.226851*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.605928*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 11.509658*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.211538*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.264045*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.394145*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.928598*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 24.077994*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 13.455648*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Ir (Multishifter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Ir()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  double omul = material->omul;
  
  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = -36.10445*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 169.0856*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 12.63649*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 168.4587*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 168.6006*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 12.62713*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 36.72283*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 36.69329*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.378321*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.423232*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 12.26825*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 12.37218*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.00601*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 12.0258*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.541984*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 165.6069*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 129.068*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 2.04229*cos(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 9.049285*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 171.4216*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 171.323*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.98783*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.72661*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.7509*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 124.7591*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 124.7676*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.529464*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.580177*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.035971*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.000332*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.279482*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.252283*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 911.255;

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(516.272*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 16.33832*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 54.29571*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 342.8432*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 342.646*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 53.92698*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 42.90644*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 43.0036*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 249.5182*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 249.5352*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 34.117856*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 34.320708*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 16.071942*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.000664*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.117928*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.009132*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 75.81894*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 336.9174*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 337.2012*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 75.76278*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 146.89132*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 146.77316*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.269926*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.539392*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 49.073*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 49.48872*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 24.01202*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 24.0516*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(38.167936*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 331.2138*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 258.136*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 8.16916*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 18.09857*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 342.646*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 35.95132*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 43.0036*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 249.5182*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 499.0704*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 17.058928*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.481062*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 32.143884*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 48.001992*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.117928*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 18.018264*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 144.4178*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 338.1712*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 25.27298*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 337.2012*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 50.50852*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 146.77316*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.539392*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 24.5365*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 74.23308*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 48.02404*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 72.1548*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-38.167936*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 331.2138*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 258.136*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 8.16916*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 36.19714*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 342.8432*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 17.97566*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 42.90644*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 499.0364*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 249.5352*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 51.176784*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 17.160354*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 48.215826*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 32.001328*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 18.235856*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.009132*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 144.4178*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 338.1712*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 50.54596*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 336.9174*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.25426*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 146.89132*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.269926*M_PI*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 73.6095*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 24.74436*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 72.03606*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 48.1032*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          if(/*(lxs+i) == 0 && */((k < int(N3/8.0)) || k > int(N3*7.0/8.0))){
            dE_core[index1] = omul*dE_core[index1];
            dE_core[index2] = omul*dE_core[index2];
            dE_core[index3] = omul*dE_core[index3];
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Ni (Multishifter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Ni()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  double omul = material->omul;
  
  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = -16.11839*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 223.4281*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 3.97842*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 223.0412*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 223.1167*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.886615*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.5071*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.51017*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.94903*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.942342*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.8734*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.894146*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.727294*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 156.1637*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 39.16039*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 2.092111*cos(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 2.812001*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 156.8333*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 156.9383*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.663776*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.992983*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.025235*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 39.01026*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 38.94427*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.727462*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.646377*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.674259*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.678564*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.193741*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.105086*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 587.2152;
          
          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(156.64156*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 16.736888*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 16.872006*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 313.6666*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 313.8766*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.982656*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 19.971932*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 20.10094*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 78.02052*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 77.88854*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.909848*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.585508*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.348518*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.357128*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.774964*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.420344*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 23.87052*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 446.0824*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 446.2334*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 23.31969*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 66.0284*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 66.04068*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 15.79612*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.769368*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.7468*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.788292*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(18.909176*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 312.3274*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 78.32078*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 8.368444*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 5.624002*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 313.8766*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.655104*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 20.10094*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 78.02052*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 155.77708*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.454924*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 15.878262*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.697036*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 16.071384*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.774964*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.840688*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 64.47356*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 446.8562*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 7.95684*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 446.2334*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 15.54646*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 66.04068*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.89806*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 23.654052*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.4936*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 23.364876*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-18.909176*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 312.3274*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 78.32078*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 8.368444*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 11.248004*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 313.6666*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.327552*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 19.971932*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 156.04104*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 77.88854*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.364772*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.292754*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.045554*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.714256*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.549928*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.420344*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 64.47356*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 446.8562*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 15.91368*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 446.0824*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.77323*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 66.0284*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 23.69418*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.884684*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 23.2404*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 15.576584*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          if(/*(lxs+i) == 0 && */((k < int(N3/8.0)) || k > int(N3*7.0/8.0))){
            dE_core[index1] = omul*dE_core[index1];
            dE_core[index2] = omul*dE_core[index2];
            dE_core[index3] = omul*dE_core[index3];
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Pd (Multishifter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Pd()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  double omul = material->omul;
  
  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = -9.18885*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 114.5172*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 2.756335*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 112.4891*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 112.4729*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.692878*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.204379*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.319446*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.52023*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.683672*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.675938*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.664205*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.689847*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 91.07588*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 26.50037*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 1.790147*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 90.20001*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 90.16154*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.743498*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.356969*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.416386*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.91209*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.90094*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.70022*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.739534*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.757776*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.855397*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 357.1646;

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(106.00148*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 10.740882*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 180.40002*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 180.32308*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.460988*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.427876*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.665544*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.82418*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 51.80188*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.80088*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.958136*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.515552*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.710794*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 16.53801*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 224.9782*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 224.9458*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 16.157268*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 36.817516*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 37.277784*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.08092*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.734688*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.351876*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.32841*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(6.759388*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 182.15176*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 53.00074*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 3.580294*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 180.32308*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.973992*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.665544*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 51.82418*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 103.60376*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.40044*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.437204*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.031104*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 11.132382*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 36.7554*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 229.0344*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 5.51267*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 224.9458*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.771512*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 37.277784*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.04046*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 16.102032*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.703752*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.98523*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-6.759388*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 182.15176*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 53.00074*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 7.160588*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 180.40002*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.486996*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.427876*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 103.64836*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 51.80188*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.20132*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 3.479068*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 10.546656*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.421588*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 36.7554*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 229.0344*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 11.02534*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 224.9782*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 5.385756*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 36.817516*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.12138*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.367344*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 16.055628*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 10.65682*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          if(/*(lxs+i) == 0 && */((k < int(N3/8.0)) || k > int(N3*7.0/8.0))){
            dE_core[index1] = omul*dE_core[index1];
            dE_core[index2] = omul*dE_core[index2];
            dE_core[index3] = omul*dE_core[index3];
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Pt (Multishifter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Pt()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  double omul = material->omul;
  
  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = -18.261*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 66.10318*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 4.290886*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 65.84877*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 65.78043*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.284049*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 17.98274*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 17.92025*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.010182*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.455881*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.371351*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.05246*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.443626*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.500909*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.017681*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.111594*sin(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.078193*sin(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.61444*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 89.68651*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 31.20643*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 1.971383*cos(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 3.566597*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 87.35806*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 87.38135*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.564238*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.164061*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.083931*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.074655*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 31.80146*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 31.84817*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 1.071898*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.37071*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.382113*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.593568*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 3.600463*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.924785*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.822556*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 379.7967;

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(124.82572*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 15.771064*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 21.399582*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 174.71612*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 174.7627*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 21.385428*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.656244*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.335724*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.44793*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 63.60292*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 63.69634*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.431388*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 13.48284*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 13.528452*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.187136*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.200926*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.69914*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.290224*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 25.745316*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 131.69754*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 131.56086*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 25.704294*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 71.93096*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 71.681*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.081456*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 17.823524*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 17.485404*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.31476*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.887252*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.001818*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.106086*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.223188*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.156386*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(6.45776*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 179.37302*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 62.41286*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 7.885532*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 7.133194*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 174.7627*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.256952*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 4.335724*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 63.60292*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 127.39268*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.431388*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.74142*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 20.292678*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.374272*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 21.602778*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.69914*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.580448*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 73.044*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 132.20636*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 8.581772*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 131.56086*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 17.136196*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 71.681*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.061092*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.911762*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 26.228106*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.10492*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 17.774504*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 27.005454*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.141448*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.669564*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.625544*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-6.45776*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 179.37302*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 62.41286*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 7.885532*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 14.266388*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 174.71612*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.128476*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.656244*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.44793*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 127.20584*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 63.69634*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 20.22426*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.764226*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 21.561408*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.401852*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 15.39828*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.290224*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 73.044*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 132.20636*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 17.163544*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 131.69754*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.568098*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 71.93096*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 2.020364*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 8.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 26.735286*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.742702*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.41968*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 26.661756*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 18.003636*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 2.035362*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 8.892752*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 6.469158*M_PI*cos(14.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          if(/*(lxs+i) == 0 && */((k < int(N3/8.0)) || k > int(N3*7.0/8.0))){
            dE_core[index1] = omul*dE_core[index1];
            dE_core[index2] = omul*dE_core[index2];
            dE_core[index3] = omul*dE_core[index3];
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
core energy Rh (Multishifter)
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_Rh()
{
  int i, j, k, isa, index, index1, index2, index3, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1, xi2, xi3;

  double b = material->b;
  double mu = material->mu;
  double omul = material->omul;
  
  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;

          fcore[index] = -23.69495*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 178.9287*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) + 7.210314*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 176.5185*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 176.5235*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.132479*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 23.90637*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 23.71257*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.24299*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.183771*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.017633*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 7.267438*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 6.729313*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 138.9095*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) - 88.26503*cos(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 1.188377*cos(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 1.357898*cos(M_PI*(16*xi1 - 8.0*xi2 - 8.0*xi3)) - 5.190815*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 141.9362*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 141.752*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.831574*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.062793*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 7.931369*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 88.16433*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 88.26113*cos(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.827993*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.490554*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.834101*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.900389*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.073984*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 1.389501*cos(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 695.5565;

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = (0.001/(mu*dslip*b))*(353.06012*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) + 9.507016*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) - 21.726368*M_PI*sin(M_PI*(16*xi1 - 8.0*xi2 - 8.0*xi3)) - 31.14489*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 283.8724*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 283.504*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 28.989444*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 32.251172*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 31.725476*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 176.32866*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 176.52226*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 19.311972*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 17.962216*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.668202*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 9.800778*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.295936*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.558004*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 43.261884*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 353.037*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 353.047*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 42.794874*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 95.62548*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 94.85028*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 28.97196*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 28.735084*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.035266*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.534876*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = (0.001/(mu*dslip*b))*(26.917252*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) - 277.819*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 176.53006*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 4.753508*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 10.863184*M_PI*sin(M_PI*(16*xi1 - 8.0*xi2 - 8.0*xi3)) + 10.38163*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 283.504*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 19.326296*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 31.725476*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 176.32866*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 353.04452*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.655986*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 26.943324*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 19.336404*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 29.402334*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 4.295936*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 11.116008*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 94.7798*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) + 357.8574*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 14.420628*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 353.047*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 28.529916*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 94.85028*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.48598*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 43.102626*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 28.070532*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 43.604628*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = (0.001/(mu*dslip*b))*(-26.917252*M_PI*sin(M_PI*(-4.0*xi2 + 4.0*xi3)) + 277.819*M_PI*sin(M_PI*(-2.0*xi2 + 2.0*xi3)) - 176.53006*M_PI*sin(M_PI*(4*xi1 - 2.0*xi2 - 2.0*xi3)) - 4.753508*M_PI*sin(M_PI*(8*xi1 - 4.0*xi2 - 4.0*xi3)) + 10.863184*M_PI*sin(M_PI*(16*xi1 - 8.0*xi2 - 8.0*xi3)) + 20.76326*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 283.8724*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 9.663148*M_PI*sin(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 32.251172*M_PI*sin(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 352.65732*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 176.52226*M_PI*sin(6.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 28.967958*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.981108*M_PI*sin(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 29.004606*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 19.601556*M_PI*sin(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 8.591872*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 5.558004*M_PI*sin(12.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 94.7798*M_PI*cos(M_PI*(-4.0*xi2 + 4.0*xi3)) - 357.8574*M_PI*cos(M_PI*(-2.0*xi2 + 2.0*xi3)) + 28.841256*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 353.037*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 14.264958*M_PI*cos(2.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 6.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 95.62548*M_PI*cos(4.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 43.45794*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) + 14.367542*M_PI*cos(8.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 4.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 42.105798*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) - 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)) - 29.069752*M_PI*cos(10.0*M_PI*(-0.5*xi2 + 0.5*xi3) + 2.0*M_PI*(xi1 - 0.5*xi2 - 0.5*xi3)));
          if(/*(lxs+i) == 0 && */((k < int(N3/8.0)) || k > int(N3*7.0/8.0))){
            dE_core[index1] = omul*dE_core[index1];
            dE_core[index2] = omul*dE_core[index2];
            dE_core[index3] = omul*dE_core[index3];
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
/* -----------------------------------------------------------------------
   core energy extended dislocations for 2 phases
   ---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::core_energy_extended_2Px()
{
  int i, j, k, isa, index, index1, index2, index3, index0, plane, num;
  int ND = dimension;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;
  double xi1,xi2,xi3;
  double c0A = material->c0;
  double c1A = material->c1;
  double c2A = material->c2;
  double c3A = material->c3;
  double c4A = material->c4;
  double a1A = material->a1;
  double a3A = material->a3;
  double b1 = material->b;

  double c0B = material2->c0;
  double c1B = material2->c1;
  double c2B = material2->c2;
  double c3B = material2->c3;
  double c4B = material2->c4;
  double a1B = material2->a1;
  double a3B = material2->a3;
  double b2 = material2->b;

  double mu1 = material->mu;
  double mu2 = material2->mu;

  c0A = c0A/(mu1*dslip*b1);
  c1A = c1A/(mu1*dslip*b1);
  c2A = c2A/(mu1*dslip*b1);
  c3A = c3A/(mu1*dslip*b1);
  c4A = c4A/(mu1*dslip*b1);
  a1A = a1A/(mu1*dslip*b1);
  a3A = a3A/(mu1*dslip*b1);

  c0B = c0B/(mu2*dslip*b2);
  c1B = c1B/(mu2*dslip*b2);
  c2B = c2B/(mu2*dslip*b2);
  c3B = c3B/(mu2*dslip*b2);
  c4B = c4B/(mu2*dslip*b2);
  a1B = a1B/(mu2*dslip*b2);
  a3B = a3B/(mu2*dslip*b2);

  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index0 = i*N2*N3 + j*N3 + k;
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
          index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

          xi1 = data_fftw[index1].re;
          xi2 = data_fftw[index2].re;
          xi3 = data_fftw[index3].re;
          if(phase[index0] == 1){
          fcore[index] = (c0B + c1B*(cos(2.0*M_PI*(xi1-xi2)) + cos(2.0*M_PI*(xi2-xi3)) + cos(2.0*M_PI*(xi3-xi1))) + c2B*(cos(2.0*M_PI*(2.0*xi1-xi2-xi3)) + cos(2.0*M_PI*(2.0*xi2-xi3-xi1)) + cos(2.0*M_PI*(2.0*xi3-xi1-xi2))) + c3B*(cos(4.0*M_PI*(xi1-xi2)) + cos(4.0*M_PI*(xi2-xi3)) + cos(4.0*M_PI*(xi3-xi1))) + c4B*(cos(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + cos(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + cos(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + cos(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + cos(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + cos(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1B*(sin(2.0*M_PI*(xi1-xi2)) + sin(2.0*M_PI*(xi2-xi3)) + sin(2.0*M_PI*(xi3-xi1))) + a3B*(sin(4.0*M_PI*(xi1-xi2)) + sin(4.0*M_PI*(xi2-xi3)) + sin(4.0*M_PI*(xi3-xi1))));

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = 2.0*M_PI*(c1B*(sin(2.0*M_PI*(xi2-xi1)) + sin(2.0*M_PI*(xi3-xi1))) + c2B*(2.0*sin(2.0*M_PI*(xi2+xi3-2.0*xi1)) + sin(2.0*M_PI*(2.0*xi2-xi3-xi1)) + sin(2.0*M_PI*(2.0*xi3-xi1-xi2))) + (2.0*c3B)*(sin(4.0*M_PI*(xi2-xi1)) + sin(4.0*M_PI*(xi3-xi1))) + c4B*(3.0*sin(2.0*M_PI*(xi2+2.0*xi3-3.0*xi1)) + 3.0*sin(2.0*M_PI*(2.0*xi2+xi3-3.0*xi1)) + 2.0*sin(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + sin(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + sin(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + 2.0*sin(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1B*(cos(2.0*M_PI*(xi1-xi2)) - cos(2.0*M_PI*(xi3-xi1))) + (2.0*a3B)*(cos(4.0*M_PI*(xi1-xi2)) - cos(4.0*M_PI*(xi3-xi1))));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = 2.0*M_PI*(c1B*(sin(2.0*M_PI*(xi1-xi2)) + sin(2.0*M_PI*(xi3-xi2))) + c2B*(sin(2.0*M_PI*(2.0*xi1-xi2-xi3)) + 2.0*sin(2.0*M_PI*(xi3+xi1-2.0*xi2)) + sin(2.0*M_PI*(2.0*xi3-xi1-xi2))) + (2.0*c3B)*(sin(4.0*M_PI*(xi1-xi2)) + sin(4.0*M_PI*(xi3-xi2))) + c4B*(sin(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + 2.0*sin(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + 3.0*sin(2.0*M_PI*(xi3+2.0*xi1-3.0*xi2)) + 3.0*sin(2.0*M_PI*(2.0*xi3+xi1-3.0*xi2)) + 2.0*sin(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + sin(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1B*(cos(2.0*M_PI*(xi2-xi3)) - cos(2.0*M_PI*(xi1-xi2))) + (2.0*a3B)*(cos(4.0*M_PI*(xi2-xi3)) - cos(4.0*M_PI*(xi1-xi2))));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = 2.0*M_PI*(c1B*(sin(2.0*M_PI*(xi2-xi3)) + sin(2.0*M_PI*(xi1-xi3))) + c2B*(sin(2.0*M_PI*(2.0*xi1-xi2-xi3)) + sin(2.0*M_PI*(2.0*xi2-xi3-xi1)) + 2.0*sin(2.0*M_PI*(xi1+xi2-2.0*xi3))) + (2.0*c3B)*(sin(4.0*M_PI*(xi2-xi3)) + sin(4.0*M_PI*(xi1-xi3))) + c4B*(2.0*sin(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + sin(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + sin(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + 2.0*sin(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + 3.0*sin(2.0*M_PI*(xi1+2.0*xi2-3.0*xi3)) + 3.0*sin(2.0*M_PI*(2.0*xi1+xi2-3.0*xi3))) + a1B*(cos(2.0*M_PI*(xi3-xi1)) - cos(2.0*M_PI*(xi2-xi3))) + (2.0*a3B)*(cos(4.0*M_PI*(xi3-xi1)) - cos(4.0*M_PI*(xi2-xi3))));
          }
          else{
          fcore[index] = (c0A + c1A*(cos(2.0*M_PI*(xi1-xi2)) + cos(2.0*M_PI*(xi2-xi3)) + cos(2.0*M_PI*(xi3-xi1))) + c2A*(cos(2.0*M_PI*(2.0*xi1-xi2-xi3)) + cos(2.0*M_PI*(2.0*xi2-xi3-xi1)) + cos(2.0*M_PI*(2.0*xi3-xi1-xi2))) + c3A*(cos(4.0*M_PI*(xi1-xi2)) + cos(4.0*M_PI*(xi2-xi3)) + cos(4.0*M_PI*(xi3-xi1))) + c4A*(cos(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + cos(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + cos(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + cos(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + cos(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + cos(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1A*(sin(2.0*M_PI*(xi1-xi2)) + sin(2.0*M_PI*(xi2-xi3)) + sin(2.0*M_PI*(xi3-xi1))) + a3A*(sin(4.0*M_PI*(xi1-xi2)) + sin(4.0*M_PI*(xi2-xi3)) + sin(4.0*M_PI*(xi3-xi1))));

          /*partial derivative wrt phase field 1*/
          dE_core[index1] = 2.0*M_PI*(c1A*(sin(2.0*M_PI*(xi2-xi1)) + sin(2.0*M_PI*(xi3-xi1))) + c2A*(2.0*sin(2.0*M_PI*(xi2+xi3-2.0*xi1)) + sin(2.0*M_PI*(2.0*xi2-xi3-xi1)) + sin(2.0*M_PI*(2.0*xi3-xi1-xi2))) + (2.0*c3A)*(sin(4.0*M_PI*(xi2-xi1)) + sin(4.0*M_PI*(xi3-xi1))) + c4A*(3.0*sin(2.0*M_PI*(xi2+2.0*xi3-3.0*xi1)) + 3.0*sin(2.0*M_PI*(2.0*xi2+xi3-3.0*xi1)) + 2.0*sin(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + sin(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + sin(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + 2.0*sin(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1A*(cos(2.0*M_PI*(xi1-xi2)) - cos(2.0*M_PI*(xi3-xi1))) + (2.0*a3A)*(cos(4.0*M_PI*(xi1-xi2)) - cos(4.0*M_PI*(xi3-xi1))));

          /*partial derivative wrt phase field 2*/
          dE_core[index2] = 2.0*M_PI*(c1A*(sin(2.0*M_PI*(xi1-xi2)) + sin(2.0*M_PI*(xi3-xi2))) + c2A*(sin(2.0*M_PI*(2.0*xi1-xi2-xi3)) + 2.0*sin(2.0*M_PI*(xi3+xi1-2.0*xi2)) + sin(2.0*M_PI*(2.0*xi3-xi1-xi2))) + (2.0*c3A)*(sin(4.0*M_PI*(xi1-xi2)) + sin(4.0*M_PI*(xi3-xi2))) + c4A*(sin(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + 2.0*sin(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + 3.0*sin(2.0*M_PI*(xi3+2.0*xi1-3.0*xi2)) + 3.0*sin(2.0*M_PI*(2.0*xi3+xi1-3.0*xi2)) + 2.0*sin(2.0*M_PI*(3.0*xi3-xi1-2.0*xi2)) + sin(2.0*M_PI*(3.0*xi3-2.0*xi1-xi2))) + a1A*(cos(2.0*M_PI*(xi2-xi3)) - cos(2.0*M_PI*(xi1-xi2))) + (2.0*a3A)*(cos(4.0*M_PI*(xi2-xi3)) - cos(4.0*M_PI*(xi1-xi2))));

          /*partial derivative wrt phase field 3*/
          dE_core[index3] = 2.0*M_PI*(c1A*(sin(2.0*M_PI*(xi2-xi3)) + sin(2.0*M_PI*(xi1-xi3))) + c2A*(sin(2.0*M_PI*(2.0*xi1-xi2-xi3)) + sin(2.0*M_PI*(2.0*xi2-xi3-xi1)) + 2.0*sin(2.0*M_PI*(xi1+xi2-2.0*xi3))) + (2.0*c3A)*(sin(4.0*M_PI*(xi2-xi3)) + sin(4.0*M_PI*(xi1-xi3))) + c4A*(2.0*sin(2.0*M_PI*(3.0*xi1-xi2-2.0*xi3)) + sin(2.0*M_PI*(3.0*xi1-2.0*xi2-xi3)) + sin(2.0*M_PI*(3.0*xi2-xi3-2.0*xi1)) + 2.0*sin(2.0*M_PI*(3.0*xi2-2.0*xi3-xi1)) + 3.0*sin(2.0*M_PI*(xi1+2.0*xi2-3.0*xi3)) + 3.0*sin(2.0*M_PI*(2.0*xi1+xi2-3.0*xi3))) + a1A*(cos(2.0*M_PI*(xi3-xi1)) - cos(2.0*M_PI*(xi2-xi3))) + (2.0*a3A)*(cos(4.0*M_PI*(xi3-xi1)) - cos(4.0*M_PI*(xi2-xi3))));
          }
        }/*ijk*/
      }
    }
  }/*plane*/
}
    /* -----------------------------------------------------------------------
    core energy sine approximation
    ---------------------------------------------------------------------*/
    /* Added by Claire 8/1/18: a 1D Sine approximations for use modeling hcp
    non-basal slip planes with GSFE that suggests partial dislocations */
    void FFTW_Slab_Conduct::core_energy_sine()
    {
      int i, j, k, isa, index, num;
      int tag;

      int ND = dimension;
      int N1 = nx;
      double size = static_cast<double>(N1);
      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NP = num_planes;
      int NS = slip_systems;
      double dslip = material->dslip;

      int NT = solve->max_iter;
      int NSI = app->stoptime;

      double c0 = material->c0;
      double c1 = material->c1;
      double c2 = material->c2;
      double c3 = material->c3;
      double c4 = material->c4;
      double a1 = material->a1;
      double a3 = material->a3;
      double isf = material->isf;
      double usf = material->usf;
      double An = material->An;
      double Cn = material->Cn;
      double b = material->b;

      double mu = material->mu;
      double ll = material->ll;
      double young = material->young;
      double xnu = material->xnu;

      c0 = c0/(mu);
      c1 = c1/(mu);
      c2 = c2/(mu);
      c3 = c3/(mu);
      c4 = c4/(mu);
      a1 = a1/(mu);
      a3 = a3/(mu);
      isf = isf/(mu*dslip*b);
      An = An/(mu*dslip*b);
      Cn = (usf - (isf/2.0))/(mu*dslip*b);
      E_core = 0;

      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          E_core += (isf*(sin(M_PI*data_fftw[index].re)*sin(M_PI*data_fftw[index].re)) + Cn*(sin(2*M_PI*data_fftw[index].re)*sin(2*M_PI*data_fftw[index].re)))/N1;
          dE_core[index] = (isf*M_PI*sin(2*M_PI*data_fftw[index].re) + Cn*2*M_PI*sin(4*M_PI*data_fftw[index].re));
        }/*ijk*/
      }/*isa*/
    }

    /* -----------------------------------------------------------------------
    core energy gamma surface
    ----------------------------------------------------------------------*/
    /* Changed by Claire 8/1/18: was core_energy_isf_usf but I
    changed to core_energy_pyrII to account for pyramidal II slip in HCP uses
    a fourier series to parameterize the 1D GSFE that governs the dissociation of
    perfect dislocations into partials with Burgers vectors in the same Direction
    with half the magnitude. The GSFE curve isn't symmetric which is why the
    parameterization is different than that for core_flag == 3, the sine approx */
    void FFTW_Slab_Conduct::core_energy_pyrII()
    {
      int i, j, k, isa, index, num;
      int tag;

      int ND = dimension;
      int N1 = nx;
      double size = static_cast<double>(N1);
      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NP = num_planes;
      int NS = slip_systems;
      double dslip = material->dslip;

      int NT = solve->max_iter;
      int NSI = app->stoptime;

      double aa0 = material->aa0;
      double aa1 = material->aa1;
      double aa2 = material->aa2;
      double aa3 = material->aa3;
      double aa4 = material->aa4;
      double bb1 = material->bb1;
      double bb2 = material->bb2;
      double bb3 = material->bb3;
      double bb4 = material->bb4;
      double isf = material->isf;
      double usf = material->usf;
      double An = material->An;
      double Cn = material->Cn;
      double b = material->b;

      double mu = material->mu;
      double ll = material->ll;
      double young = material->young;
      double xnu = material->xnu;

      aa0 = aa0/(mu*dslip*b);
      aa1 = aa1/(mu*dslip*b);
      aa2 = aa2/(mu*dslip*b);
      aa3 = aa3/(mu*dslip*b);
      aa4 = aa4/(mu*dslip*b);
      bb1 = bb1/(mu*dslip*b);
      bb2 = bb2/(mu*dslip*b);
      bb3 = bb3/(mu*dslip*b);
      bb4 = bb4/(mu*dslip*b);
      isf = isf/(mu*dslip*b);
      An = An/(mu*dslip*b);
      Cn = (usf - (isf/2.0))/(mu*dslip*b);

      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          fcore[index] = (aa0 +aa1*cos(2.0*M_PI*data_fftw[index].re)+bb1*sin(2.0*M_PI*data_fftw[index].re)+aa2*cos(4.0*M_PI*data_fftw[index].re)+bb2*sin(4.0*M_PI*data_fftw[index].re)+aa3*cos(6.0*M_PI*data_fftw[index].re)+bb3*sin(6.0*M_PI*data_fftw[index].re)+aa4*cos(8.0*M_PI*data_fftw[index].re)+bb4*sin(8.0*M_PI*data_fftw[index].re));

          MPI_Reduce(&fcore[index], &f_core[index], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }/*ijk*/

        dE_core[index] = ((2.0*M_PI)*(bb1*cos(2.0*M_PI*data_fftw[index].re)-aa1*sin(2.0*M_PI*data_fftw[index].re))+(4.0*M_PI)*(bb2*cos(4.0*M_PI*data_fftw[index].re)-aa2*sin(4.0*M_PI*data_fftw[index].re))+(6.0*M_PI)*(bb3*cos(6.0*M_PI*data_fftw[index].re)-aa3*sin(6.0*M_PI*data_fftw[index].re))+(8.0*M_PI)*(bb4*cos(8.0*M_PI*data_fftw[index].re)-aa4*sin(8.0*M_PI*data_fftw[index].re)));
      }/*isa*/
    }

    /*----------------------------------------
    The stress-state dependet USFE magnitude: USFE(MRSSP angle X, shear stress tau)
    -----------------------------------------*/
    /* added by Hyojung 12/2020*/

    void FFTW_Slab_Conduct::core_energy_USFE_angle_tau()
    {
      int i, j, k, isa, index, num;
      int tag;
      int N1 = nx;
      double size = static_cast<double>(N1);
      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NS = slip_systems;
      double dslip = material->dslip;
      double usf = material->usf;
      double An ;
      double b = material->b;
      double mu = material->mu;

      double a_slope = material->a_slope;
      double a_b = material->a_b;
      double b_b = material->b_b;
      double c_slope= material->c_slope;
      double c_b = material->c_b;
      double angle_to_110 = app->angle_to_110;

      double gsfe_max = (a_slope*sigma_rot[2][1]*mu*1.0E-9+a_b)*sin(b_b*angle_to_110*M_PI/180.)+(c_slope*sigma_rot[2][1]*mu*1.0E-9+c_b) ;
      //  An = usf/(mu*dslip*b);
      gsfe_max = gsfe_max*1.0E-3;  // unit: mJ/m^2 to J/m^
      An = gsfe_max/(mu*dslip*b);

      if(me==0){
        if (logfile) fprintf(logfile,"core flag=6, usfe(angle,tau): %lf\n", gsfe_max);
        if (logfile) fprintf(logfile,"sigma: %lf, angle_to_110: %lf\n", sigma_rot[1][2], angle_to_110);
      }

      E_core = 0;

      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          E_core += An*(sin(M_PI*data_fftw[index].re)*sin(M_PI*data_fftw[index].re))/N1;   //To make it general fft->data_fftw has to be general
          dE_core[index] = An*M_PI*sin(2.0*M_PI*data_fftw[index].re);
        }/*ijk*/
      }/*isa*/
    }

    /* -----------------------------------------------------------------------
    Project core energy for output
    ---------------------------------------------------------------------*/
    void FFTW_Slab_Conduct::project_core_energy()
    {
      int index, index1, index2, index3, indexdx;
      int counter, marker, indexmin, tag, countSF, cSF, count;
      double dx, mpidel, p, totAR, sfAR;
      MPI_Status status;
      int ND = dimension;
      int N1 = nx;
      double size = static_cast<double>(nx);
      int lN1 = local_x;
      int lxs = local_x_start;
      int N2 = local_y;
      int N3 = local_z;
      int NP = num_planes;
      int NS = slip_systems;
      double dslip = material->dslip;
      double b = material->b;

      dx = size/N1;
      tag = 1;
      p = M_PI/(sqrt(3.0)*(b/1.0E-10));

      for(int plane=0;plane<NP;plane++){
        for(int i=0;i<lN1;i++)
        for(int j=0;j<N2;j++)
        for(int k=0;k<N3;k++){
          if(NS == 1 || NS == 2){
            index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
            delta[index] = data_fftw[index].re;
          } // NS=1 || 2
          else{
            index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
            //indexdx = ((lxs+i)+1)*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
            index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
            index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
            index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

            delta[index] = (data_fftw[index1].re*xb[0][0] + data_fftw[index2].re*xb[1][0]
              + data_fftw[index3].re*xb[2][0])*xb[1][0]
              + (data_fftw[index1].re*xb[0][1] + data_fftw[index2].re*xb[1][1]
                + data_fftw[index3].re*xb[2][1])*xb[1][1]
                + (data_fftw[index1].re*xb[0][2] + data_fftw[index2].re*xb[1][2] + data_fftw[index3].re*xb[2][2])*xb[1][2];

              }
            } //ijk
          }//plane

          for(int plane=0; plane<NP; plane++){
            for(int i=0;i<lN1;i++){
              int j = N2/2;
              int k = N3/2;
              index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
              indexdx = (i+1)*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
              if((i+1) == lN1){
                indexmin = 0*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
                if(me != 0){
                  //	  fft_send(&delta[indexmin], 1, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);   // This should be a send function in the FFT class particular to the
                  MPI_Send(&delta[indexmin], 1, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
                }
                if((lxs+(i+1)) == N1){
                  mpidel = delta[index];
                  ddelta[index] = (mpidel - delta[index])/dx;
                  break;
                }
                //	fft_recv(&mpidel, 1, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, status);
                MPI_Recv(&mpidel, 1, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
                ddelta[index] = (mpidel - delta[index])/dx;
              }
              else{
                ddelta[index] = (delta[indexdx] - delta[index])/dx;
              }
            }
          }

        }
 
  /* -----------------------------------------------------------------------
    core energy for perfect dislocations
    ---------------------------------------------------------------------*/
  void FFTW_Slab_Conduct::gb_energy_proportional()
{
  int i, j, k, isa, index, num, index0;
  int tag;

  int ND = dimension;
  int N1 = nx;
  double size = static_cast<double>(N1);
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;

  int NT = solve->max_iter;
  int NSI = app->stoptime;

  double c0 = material->c0;
  double c1 = material->c1;
  double c2 = material->c2;
  double c3 = material->c3;
  double c4 = material->c4;
  double a1 = material->a1;
  double a3 = material->a3;
  double isf = material->isf;
  double usf = material->usf;
  double An = material->An;
  double Cn = material->Cn;
  double b = material->b;

  double mu = material->mu;
  double ll = material->ll;
  double young = material->young;
  double xnu = material->xnu;

  if(app->nphase > 1){
  double muinc = material2->mu;
  double pbef1  = material->pbef;
  double pbef2 = material2->pbef;

  double gbe = material2->gbe;
  double gbef = material2->gbef;

  c0 = c0/(mu);
  c1 = c1/(mu);
  c2 = c2/(mu);
  c3 = c3/(mu);
  c4 = c4/(mu);
  a1 = a1/(mu);
  a3 = a3/(mu);
  isf = isf/(mu*dslip*b);
  An = An/(mu*dslip*b);
  Cn = (usf - (isf/2.0))/(mu*dslip*b);

  E_core = 0;
  // temp 2phase core energy
  double An1;
  
  An = pbef1*usf/(mu*dslip*b);
  // An1 = usf/(muinc*dslip*b);
  An1 = pbef2*usf/(mu*dslip*b);

  // if(me == 0){
  //   if (screen) fprintf(screen, "pbef1 = %lf \n",pbef1);
  //   if (screen) fprintf(screen, "pbef2 = %lf \n",pbef2);
  //   if (screen) fprintf(screen, "An1 = %lf \n", An);
  //   if (screen) fprintf(screen, "An2 = %lf\n", An1);
  // }

  for(isa=0;isa<NS;isa++){
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
      dE_gb[index] = 0.0; // no gb energy inside the grain
      if (app->nphase>1) {
        index0 = i*N2*N3 + j*N3 + k;
        if (phase[index0] == 1) {
          // dE_gb[index] = dE_core[index]*gbe/(muinc*dslip*b*An1);
          dE_gb[index] = dE_core[index]*(gbe - usf)/(muinc*dslip*b*An1)*gbef;
          // dE_gb[index] = dE_core[index]*gbe/(muinc*dslip*b); // factor
          // if (fft->me == 0) {
          //   if (screen) fprintf(screen,"dE_core[%d] = %lf,  dE_gb[%d] = %lf\n", index, dE_core[index], index, dE_gb[index]);
          //   if (logfile) fprintf(logfile,"dE_core[%d] = %lf,  dE_gb[%d] = %lf\n", index, dE_core[index], index, dE_gb[index]);
          // }
        } 
      }

    }/*ijk*/
  }/*isa*/
  }
}

  /* -----------------------------------------------------------------------
  grain boundary energy extended 
  ---------------------------------------------------------------------*/
  /* This is just the extended core energy with a factor for the 
  grain boundary energy */
  void FFTW_Slab_Conduct::gb_energy_extended()
  {
    int i, j, k, isa, index, index0, index1, index2, index3, plane, num;
    int tag;
    //  int counter, marker, tag, countSF, cSF, count;

    int ND = dimension;
    int N1 = nx;
    double size = static_cast<double>(N1);
    int lN1 = local_x;
    int lxs = local_x_start;
    int N2 = local_y;
    int N3 = local_z;
    int NP = num_planes;
    int NS = slip_systems;
    double dslip = material->dslip;

    int NT = solve->max_iter;
    int NSI = app->stoptime;

    double c0 = material->c0;
    double c1 = material->c1;
    double c2 = material->c2;
    double c3 = material->c3;
    double c4 = material->c4;
    double a1 = material->a1;
    double a3 = material->a3;
    double isf = material->isf;
    double usf = material->usf;
    double An = material->An;
    double Cn = material->Cn;
    double b = material->b;

    double mu = material->mu;
    double ll = material->ll;
    double young = material->young;
    double xnu = material->xnu;

    double muinc = material2->mu;
    double pbef1  = material->pbef;
    double pbef2 = material2->pbef;

    double gbe = material2->gbe;

    double An1;

    c0 = c0/(mu);
    c1 = c1/(mu);
    c2 = c2/(mu);
    c3 = c3/(mu);
    c4 = c4/(mu);
    a1 = a1/(mu);
    a3 = a3/(mu);
    isf = isf/(mu*dslip*b);
    Cn = (usf - (isf/2.0))/(mu*dslip*b);

    An = pbef1*usf/(mu*dslip*b);
    An1 = pbef2*usf/(mu*dslip*b);

    /*as of now the planes need to be in the correct order initially for this to work
    also make sure to call this subroutine before the FFTW functions are called so that data_fftw = xi[na0] or double check that this is true.
    subroutine calculates core energy for each time step, need to called every time step and initialize every time step.*/

    for(isa=0;isa<NS;isa++){
      for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
      for(k=0;k<N3;k++){
        index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
        data_core[index].re = data_fftw[index].re;
        data_core[index].im = data_fftw[index].im;
      }/*ijk*/
    }/*isa*/

    for(plane=0;plane<NP;plane++){
      for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
      for(k=0;k<N3;k++){
        index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
        index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
        index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
        index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

        /*no derivatives taken use to calculate E_core*/

        fcore[index] = (c0 + c1*(cos(2.0*M_PI*(data_core[index1].re-data_core[index2].re)) + cos(2.0*M_PI*(data_core[index2].re-data_core[index3].re)) + cos(2.0*M_PI*(data_core[index3].re-data_core[index1].re))) + c2*(cos(2.0*M_PI*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + cos(2.0*M_PI*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + cos(2.0*M_PI*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + c3*(cos(4.0*M_PI*(data_core[index1].re-data_core[index2].re)) + cos(4.0*M_PI*(data_core[index2].re-data_core[index3].re)) + cos(4.0*M_PI*(data_core[index3].re-data_core[index1].re))) + c4*(cos(2.0*M_PI*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + cos(2.0*M_PI*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + cos(2*M_PI*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + cos(2.0*M_PI*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + cos(2.0*M_PI*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + cos(2.0*M_PI*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(sin(2.0*M_PI*(data_core[index1].re-data_core[index2].re)) + sin(2.0*M_PI*(data_core[index2].re-data_core[index3].re)) + sin(2.0*M_PI*(data_core[index3].re-data_core[index1].re))) + a3*(sin(4.0*M_PI*(data_core[index1].re-data_core[index2].re)) + sin(4.0*M_PI*(data_core[index2].re-data_core[index3].re)) + sin(4.0*M_PI*(data_core[index3].re-data_core[index1].re))))/(dslip*b);

        MPI_Reduce(&fcore[index], &f_core[index], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        //            E_core += (f_core[index]/nsize);

        /*partial derivative wrt phase field 1*/
        df1core[index] = ((2.0*M_PI)/(dslip*b))*(c1*(sin(2.0*M_PI*(data_core[index2].re-data_core[index1].re)) + sin(2.0*M_PI*(data_core[index3].re-data_core[index1].re))) + c2*(2.0*sin(2.0*M_PI*(data_core[index2].re+data_core[index3].re-2.0*data_core[index1].re)) + sin(2.0*M_PI*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + sin(2.0*M_PI*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + (2.0*c3)*(sin(4.0*M_PI*(data_core[index2].re-data_core[index1].re)) + sin(4.0*M_PI*(data_core[index3].re-data_core[index1].re))) + c4*(3.0*sin(2.0*M_PI*(data_core[index2].re+2.0*data_core[index3].re-3.0*data_core[index1].re)) + 3.0*sin(2.0*M_PI*(2.0*data_core[index2].re+data_core[index3].re-3.0*data_core[index1].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + sin(2.0*M_PI*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + sin(2.0*M_PI*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(cos(2.0*M_PI*(data_core[index1].re-data_core[index2].re)) - cos(2.0*M_PI*(data_core[index3].re-data_core[index1].re))) + (2.0*a3)*(cos(4.0*M_PI*(data_core[index1].re-data_core[index2].re)) - cos(4.0*M_PI*(data_core[index3].re-data_core[index1].re))));

        /*partial derivative wrt phase field 2*/
        df2core[index] = ((2.0*M_PI)/(dslip*b))*(c1*(sin(2.0*M_PI*(data_core[index1].re-data_core[index2].re)) + sin(2.0*M_PI*(data_core[index3].re-data_core[index2].re))) + c2*(sin(2.0*M_PI*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + 2.0*sin(2.0*M_PI*(data_core[index3].re+data_core[index1].re-2.0*data_core[index2].re)) + sin(2.0*M_PI*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + (2.0*c3)*(sin(4.0*M_PI*(data_core[index1].re-data_core[index2].re)) + sin(4.0*M_PI*(data_core[index3].re-data_core[index2].re))) + c4*(sin(2.0*M_PI*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + 3.0*sin(2.0*M_PI*(data_core[index3].re+2.0*data_core[index1].re-3.0*data_core[index2].re)) + 3.0*sin(2.0*M_PI*(2.0*data_core[index3].re+data_core[index1].re-3.0*data_core[index2].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + sin(2.0*M_PI*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(cos(2.0*M_PI*(data_core[index2].re-data_core[index3].re)) - cos(2.0*M_PI*(data_core[index1].re-data_core[index2].re))) + (2.0*a3)*(cos(4.0*M_PI*(data_core[index2].re-data_core[index3].re)) - cos(4.0*M_PI*(data_core[index1].re-data_core[index2].re))));

        /*partial derivative wrt phase field 3*/
        df3core[index] = ((2.0*M_PI)/(dslip*b))*(c1*(sin(2.0*M_PI*(data_core[index2].re-data_core[index3].re)) + sin(2.0*M_PI*(data_core[index1].re-data_core[index3].re))) + c2*(sin(2.0*M_PI*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + sin(2.0*M_PI*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + 2.0*sin(2.0*M_PI*(data_core[index1].re+data_core[index2].re-2.0*data_core[index3].re))) + (2.0*c3)*(sin(4.0*M_PI*(data_core[index2].re-data_core[index3].re)) + sin(4.0*M_PI*(data_core[index1].re-data_core[index3].re))) + c4*(2.0*sin(2.0*M_PI*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + sin(2.0*M_PI*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + sin(2.0*M_PI*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + 2.0*sin(2.0*M_PI*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + 3.0*sin(2.0*M_PI*(data_core[index1].re+2.0*data_core[index2].re-3.0*data_core[index3].re)) + 3.0*sin(2.0*M_PI*(2.0*data_core[index1].re+data_core[index2].re-3.0*data_core[index3].re))) + a1*(cos(2.0*M_PI*(data_core[index3].re-data_core[index1].re)) - cos(2.0*M_PI*(data_core[index2].re-data_core[index3].re))) + (2.0*a3)*(cos(4.0*M_PI*(data_core[index3].re-data_core[index1].re)) - cos(4.0*M_PI*(data_core[index2].re-data_core[index3].re))));
      }/*ijk*/
    }/*plane*/

    for(plane=0;plane<NP;plane++){
      for(i=0;i<lN1;i++)
      for(j=0;j<N2;j++)
      for(k=0;k<N3;k++){
        index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
        index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
        index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
        index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

        dE_gb[index1] = 0.0;
        dE_gb[index2] = 0.0;
        dE_gb[index3] = 0.0;
        if (app->nphase>1) {
          index0 = i*N2*N3 + j*N3 + k;
          if (phase[index0] == 1) {
            dE_gb[index1] = df1core[index]*gbe/(muinc*dslip*b*An1);
            dE_gb[index2] = df1core[index]*gbe/(muinc*dslip*b*An1);
            dE_gb[index3] = df1core[index]*gbe/(muinc*dslip*b*An1);
            // if (fft->me == 0) {
            //   if (screen) fprintf(screen,"dE_core[%d] = %lf,  dE_gb[%d] = %lf\n", index, dE_core[index], index, dE_gb[index]);
            //   if (logfile) fprintf(logfile,"dE_core[%d] = %lf,  dE_gb[%d] = %lf\n", index, dE_core[index], index, dE_gb[index]);
            // }
          }
        }
      }/*ijk*/
    }/*plane*/
    MPI_Barrier(MPI_COMM_WORLD);

  }
 
/* -----------------------------------------------------------------------
no grain boundary lattice energy 
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::no_gb_energy()
{
  int i, j, k, isa, index, index0, index1, index2, index3, plane, num;
  int tag;

  int ND = dimension;
  int N1 = nx;
  double size = static_cast<double>(N1);
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  double dslip = material->dslip;

  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
      index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + plane*lN1*N2*N3*3;
      index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + plane*lN1*N2*N3*3;
      index3 = i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + plane*lN1*N2*N3*3;

      dE_gb[index1] = 0.0;
      dE_gb[index2] = 0.0;
      dE_gb[index3] = 0.0;
    }/*ijk*/
  }/*plane*/
  MPI_Barrier(MPI_COMM_WORLD);

}

        /* -----------------------------------------------------------------------
        Initialize inner loop
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::init_loop()
        {
          int lN1 = local_x;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;
          int NP = num_planes;
          int ND = dimension;

          for(int i=0;i<lN1*N2*N3*2;i++){
            xi_sum[0][i] = 0.0;
          }

          E_core = 0.0;  //Core Energy for each time step.

          for(int i=0;i<lN1*N2*N3*NS;i++){
            dE_core[i] = 0.0;
          }

          for(int i=0;i<lN1*N2*N3*NP;i++){
            fcore[i] = 0.0;
            df1core[i] = 0.0;
            df2core[i] = 0.0;
            df3core[i] = 0.0;
            delta[i] = 0.0;
            ddelta[i] = 0.0;
          }

          // for(int i=0;i<ND;i++)
          //   for(int j=0;j<ND;j++){
          //     local_sigma[i][j] = 0.0;
          //     avepst[i][j] = 0.0;
          //     avepsts[i][j] = 0.0;
          //     avepsd[i][j] = 0.0;
          //     aveps[i][j] = 0.0;
          //     ave_eps[i][j] = 0.0;
          //     ave_sigma[i][j] = 0.0;
          //   }
        }
        /* -----------------------------------------------------------------------
        Compute Order parameter that minimizes internal energy
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::internal_energy()
        {
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;
          int ND = dimension;
          int index=0, index2=0, na=0, nb=0, nb0=0, nc=0;
          int is, isa, isb, ida, idb, i, j, k, l, m, n, u, v, na0, na02, na1, nb1, psys, ptindex;
          double mu = material->mu;
          double b = material->b;
          double dslip = material->dslip;

          for (i=0; i<lN1*N2*N3*NS; i++){
            data_real[i].re = 0;
            data_real[i].im = 0;
          }

          for(isa=0;isa<NS;isa++){
            for(isb=0;isb<NS;isb++){
              for(i=0;i<lN1;i++)
              for(j=0;j<N2;j++)
              for(k=0;k<N3;k++){
                index  = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
                index2 = i*N2*N3 + j*N3 + k + isb*lN1*N2*N3;
                nb     = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3 + isb*lN1*N2*N3*NS;
                data_real[index].re += data_fftw[index2].re * BB[nb];
                data_real[index].im += data_fftw[index2].im * BB[nb];
              }
            }
          }

          if (app->nphase>1) {
            // for(ida=0;ida<ND;ida++){
            //   for(idb=0;idb<ND;idb++){
            //     for(i=0;i<lN1;i++) {
            //       for(j=0;j<N2;j++) {
            //         for(k=0;k<N3;k++) {
            //           na0 = 2*(i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
            //           na1 = na0+1;
            //           index = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
            //           vstrain_complex[index].re = vstrain[na0];
            //           vstrain_complex[index].im = vstrain[na1];
            //         }
            //       }
            //     }
            //   }
            // }
            for(ida=0;ida<ND*ND;ida++){
              psys = ida*lN1*N2*N3;
              fftwnd_mpi(plan, 1, vstrain_complex+psys, work_strain, FFTW_NORMAL_ORDER);
            }
            for(i=0;i<lN1;i++){ // Multiply virtual strain by DD
              for(j=0;j<N2;j++){
                for(k=0;k<N3;k++) {
                  for(isa=0;isa<NS;isa++){
                    for(ida=0;ida<ND;ida++){
                      for(idb=0;idb<ND;idb++){
                        index  = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
                        index2 = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
                        nb = i*N2*N3 + j*N3 + k+ isa*lN1*N2*N3 + ida*lN1*N2*N3*NS + idb*lN1*N2*N3*NS*ND;
                        data_real[index].re += vstrain_complex[index2].re * DD[nb];
                        data_real[index].im += vstrain_complex[index2].im * DD[nb];
                      }
                    }
                  }
                }
              }
            }
          }
          //Compute divergence of electric current
          FFTW_Slab_Conduct::divcurrent();

        }
/* -----------------------------------------------------------------------
   Compute VS that minimizes internal energy
   ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::internal_energy_vs()
        {if (app->nphase>1) {
          int lN1 = local_x;
          int lxs = local_x_start;
          int N1 = nx;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;
          int ND = dimension;
          int nsize = N1*N2*N3;
          int index=0, index2=0, nb=0;
          int isa, isb, ida, idb, i, j, k, na0, na1, psys;
          int ia, ib, ic, id;
          int kpa;
          for (i=0; i<lN1*N2*N3*ND*ND; i++){
            elastic_energy_den_vs[i].re = 0;
            elastic_energy_den_vs[i].im = 0;
          }
          for(i=0;i<lN1;i++){ // Multiply order parameter by DD
            for(j=0;j<N2;j++){
              for(k=0;k<N3;k++) {
                for(isa=0;isa<NS;isa++){
                  for(ida=0;ida<ND;ida++){
                    for(idb=0;idb<ND;idb++){
                      index  = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
                      nb = i*N2*N3 + j*N3 + k+ isa*lN1*N2*N3 + ida*lN1*N2*N3*NS + idb*lN1*N2*N3*NS*ND;
                      index2 = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
                      elastic_energy_den_vs[index].re += data_fftw[index2].re * DD[nb];
                      elastic_energy_den_vs[index].im += data_fftw[index2].im * DD[nb];
                    }
                  }
                }
              }
            }
          }
          for(i=0;i<lN1;i++){ // Multiply virtual strain by AA
            for(j=0;j<N2;j++){
              for(k=0;k<N3;k++) {
                for(ia=0;ia<ND;ia++){
                  for(ib=0;ib<ND;ib++){
                    for(ic=0;ic<ND;ic++){
                      for(id=0;id<ND;id++){
                        nb = i*N2*N3 + j*N3 + k + ia*lN1*N2*N3 + ib*lN1*N2*N3*ND + ic*lN1*N2*N3*ND*ND + id*lN1*N2*N3*ND*ND*ND;
                        index  = i*N2*N3 + j*N3 + k + ia*lN1*N2*N3 + ib*lN1*N2*N3*ND;
                        index2  = i*N2*N3 + j*N3 + k + ic*lN1*N2*N3 + id*lN1*N2*N3*ND;
                        elastic_energy_den_vs[index].re += vstrain_complex[index2].re * AA[nb];
                        elastic_energy_den_vs[index].im += vstrain_complex[index2].im * AA[nb];
                      }
                    }
                  }
                }
              }
            }
          }
          for(ida=0;ida<ND*ND;ida++){
            psys = ida*lN1*N2*N3;
            fftwnd_mpi(iplan, 1, elastic_energy_den_vs+psys, work, FFTW_NORMAL_ORDER);
          }
          for (i=0; i<lN1*N2*N3*ND*ND; i++){
            elastic_energy_den_vs[i].re = elastic_energy_den_vs[i].re / ((double) nsize);
            elastic_energy_den_vs[i].im = elastic_energy_den_vs[i].im / ((double) nsize);
          }
          for(i=0;i<lN1;i++){ // Multiply virtual strain by MM
            for(j=0;j<N2;j++){
              for(k=0;k<N3;k++) {
                for(ia=0;ia<ND;ia++){
                  for(ib=0;ib<ND;ib++){
                    for(ic=0;ic<ND;ic++){
                      for(id=0;id<ND;id++){
                        nb = i*N2*N3 + j*N3 + k ;
                        index  = i*N2*N3 + j*N3 + k + ia*lN1*N2*N3 + ib*lN1*N2*N3*ND;
                        index2 = 2*(i*N2*N3 + j*N3 + k + ic*lN1*N2*N3 + id*lN1*N2*N3*ND);
                        kpa = ((int)phase[nb]) - 1;

                        if (kpa%1!=0){
  if(me == 0){
    if (logfile) fprintf(logfile,"Not integer kpa = %lf\n", kpa);
    if (screen) fprintf(screen,"Not integer kpa = %lf\n", kpa);
  }
}
                        if (kpa>=0) {
                          elastic_energy_den_vs[index].re += vstrain[index2] * MM[kpa][ia][ib][ic][id];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        }

		/* -----------------------------------------------------------------------
		   Compute Fourier space pieces of divergence of electric current from according order parameter
		   ---------------------------------------------------------------------*/
		void FFTW_Slab_Conduct::divcurrent()
		{
		  int lN1 = local_x;
		  int N2 = local_y;
		  int N3 = local_z;
		  int index=0;
		  int nsize = lN1*N2*N3;
		  
		  double fk[dimension];
		  double fk2;
		  
		  //~ double epsilonzero = material->epsilonzero;
		  double ElEx[3];
		  //~ double condsigma = material->condsigma;
		  //~ double condsigma_array[lN1];
		  
		  for (int l=0;l<dimension;l++) ElEx[l] = material->El[l];
		  
		  for (int i=0; i<lN1*N2*N3; i++){
			conddata_bpart[i].re = 0;
			conddata_bpart[i].im = 0;
			jcurrent_fftw[i].re = 0; // omit external part condsigma*ElEx here; // TODO: make this a vector
			jcurrent_fftw[i].im = 0;
		  }
		  
		  for(int i=0;i<lN1;i++)
			for(int j=0;j<N2;j++)
			  for(int k=0;k<N3;k++){
				index  = i*N2*N3 + j*N3 + k;
				fk[0] = fx[index];
				fk[1] = fy[index];
				fk[2] = fz[index];
				fk2 = fk[0]*fk[0]+fk[1]*fk[1]+fk[2]*fk[2];
				if (fk2>0.)
				{
				// calculate for y-direction for now; TODO: generalize by turning scalars into vectors and tensors
				jcurrent_fftw[index].re -= conddata_fftw[index].im*fk[1]/fk2; // tilde{a}_j in working notes
				jcurrent_fftw[index].im += conddata_fftw[index].re*fk[1]/fk2;
				conddata_bpart[index].re += conddata_fftw[index].re*fk[1]*fk[1]/fk2; // tilde{b}_{ij} in working notes
				conddata_bpart[index].im += conddata_fftw[index].im*fk[1]*fk[1]/fk2; // TODO: include non-x components
			  }}
		}
/* -----------------------------------------------------------------------
   Compute total energy
   ---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::energy_calc()
{
  int lN1 = local_x;
  int lxs = local_x_start;
  int N1 = nx;
  int N2 = local_y;
  int N3 = local_z;
  int NS = slip_systems;
  int NP = num_planes;
  int ND = dimension;
  int index=0, index2=0, na=0, nb=0, nb0=0, nc=0;
  int is, isa, isb, ida, idb, u, v, na0, na02, na1, nb1, psys, plane;
  double mu = material->mu;
  double b = material->b;
  double dslip = material->dslip;
  int i, j, k, l,m,n,ptindex,index1,index3,index4,index5,index6;
  double int_ave=0.0, ext_ave=0.0, lat_ave=0.0;
  int fcount = 0;

  for (i=0; i<lN1*N2*N3*NS; i++){
    data_real[i].re = 0;
    data_real[i].im = 0;
  }

  for (i=0; i<lN1*N2*N3; i++){
    data_real_2[i].re = 0;
    data_real_2[i].im = 0;
    tau_2[i] = 0;

    E_int[i] = 0;
    E_ext[i] = 0;
    E_lat[i] = 0;
  }

  for(is=0;is<NS;is++){
    for (ida=0; ida<ND; ida++){
      for (idb=0; idb<ND; idb++){
  for(i=0;i<lN1;i++){
    for(j=0;j<N2;j++){
      for(k=0;k<N3;k++){
        nb = 2*(k + j*N3 + i*N3*N2+ is*lN1*N2*N3);
        ptindex = k + j*N3 + i*N3*N2;
        index = k + j*N3 + i*N3*N2 + is*lN1*N2*N3;
        tau_2[ptindex] += eps[is][ida][idb] * xi[0][nb] * sigma_rot[ida][idb]*mu*b*dslip;
      }
    }
  }
      }
    }
  }

  for (ida=0;ida<ND;ida++){
    for (idb=0;idb<ND;idb++){
      for(i=0;i<lN1;i++){
        for(j=0;j<N2;j++){
          for(k=0;k<N3;k++){
            na = 2*(i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
            nb = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
            nb0 = i*N2*N3 + j*N3 + k;
            for (m=0;m<ND;m++){
              for (l=0;l<ND;l++){
                na0 = 2*(i*N2*N3 + j*N3 + k + m*lN1*N2*N3 + l*lN1*N2*N3*ND);
                index = i*N2*N3 + j*N3 + k + m*lN1*N2*N3 + l*lN1*N2*N3*ND;
                data_real_2[nb0].re += (0.5*b*dslip)*C[ida][idb][m][l]*(data_strain[index].re-data_epsd[na0]-avepsts[m][l])*(data_strain[nb].re-data_epsd[na]-avepsts[ida][idb]);
                }
              }
            }
          }
        }
      }
    }

  for(plane=0;plane<NP;plane++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + plane*lN1*N2*N3;
          ptindex = i*N2*N3 + j*N3 + k;

          E_lat[ptindex] += fcore[index];
        }
      }
    }
  }


  for(i=0;i<lN1;i++){
    for(j=0;j<N2;j++){
      for(k=0;k<N3;k++){
        ptindex = i*N2*N3 + j*N3 + k;
        index1 = 2*(i*N2*N3 + j*N3 + k + 0*lN1*N2*N3);
        index2 = 2*(i*N2*N3 + j*N3 + k + 1*lN1*N2*N3);

        E_int[ptindex] = data_real_2[ptindex].re;
        E_ext[ptindex] = tau_2[ptindex];

        int_ave += E_int[ptindex];
        ext_ave += E_ext[ptindex];
        lat_ave += E_lat[ptindex];
    }
  }
}
  MPI_Allreduce(&int_ave, &intave, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ext_ave, &extave, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&lat_ave, &latave, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // if(ecount != 0){
  //   intave /= ecount;
  //   extave /= ecount;
  //   latave /= ecount;
  // }
  // if(me == 0) printf("ecount %d\n",ecount);
}

/* -----------------------------------------------------------------------
Update order parameter
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::update_order_parameter()
{
  int lN1 = local_x;
  int lxs = local_x_start;
  int N1 = nx;
  int N2 = local_y;
  int N3 = local_z;
  int NS = slip_systems;
  int NP = num_planes;
  int index=0, index1=0, index2=0, index3=0, na=0, na0=0, na1=0, ptindex = 0, plnindex = 0, indexp = 0;
  double xinormlocal=0.0;           // stores the norm of the increment of xi
  double xirep=0.0, xiimp=0.0;      // store previous re and im values of xi
  double xinormp=xinorm;
  double xi_ave=0.0;                // Local average of the order parameter
  int nsize = N1*N2*N3;
  int i,j,k,isa,ip;
  double mu = material->mu;
  double b = material->b;
  double dslip = material->dslip;

  for(isa=0;isa<NS;isa++){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
	for(k=0;k<N3;k++){
          na0 = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          index1 = i*N2*N3 + j*N3 + k + 0*lN1*N2*N3;
          index2 = i*N2*N3 + j*N3 + k + 1*lN1*N2*N3;
          ptindex = i*N2*N3 + j*N3 + k;
          plnindex = i*N2*N3 + j*N3 + k + (isa/3)*lN1*N2*N3;
          na = 2*(i*N2*N3 + j*N3 + k);
          na1 = na0+1;

          //Ginzburg-Landau Equation for real and imag parts
            xirep = xi[0][na0];
            xiimp = xi[0][na1];
          if(xo[index] <= 0.0){ // kMC solver has negative xo
            xi[0][na0] = xi[0][na0]-((app->CD*app->timestep)*(data_real[index].re/(nsize) - tau[isa] + dE_core[index]));
            xi[0][na1] = xi[0][na1]-((app->CD*app->timestep)*(data_real[index].im/(nsize)));
            xi_sum[0][na] += xi[0][na0];
            xi_sum[0][na+1] += xi[0][na1];
          }
          xinormlocal += (xi[0][na0] - xirep)*(xi[0][na0] - xirep) +
          (xi[0][na1] - xiimp)*(xi[0][na1] - xiimp);
          // data_fftw[index].re = xi[0][na0];
          // data_fftw[index].im = xi[0][na1];
          xi_ave += xi[0][na0];
        }
      }
    }
  }
if(app->dflag == 1){
for(ip=0;ip<NP;ip++){
for(i=0;i<lN1;i++){
for(j=0;j<N2;j++){
    for(k=0;k<N3;k++){
      index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
      index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
      index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
      na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
      na1=na0+1;
      dreg[0][na0] = 2.0*xi[0][index1] - xi[0][index2] - xi[0][index3];
      dreg[0][na1] = 3.0*xi[0][index2] - 3.0*xi[0][index3];
      dregp[0][na0] = dreg[0][na0]*cos(M_PI/6.0)+dreg[0][na1]*sin(M_PI/6.0); // Trailing partial
      dregp[0][na1] = dreg[0][na0]*cos(M_PI/6.0)-dreg[0][na1]*sin(M_PI/6.0); // Leading partial        
      }
    }
  }
}
}
  MPI_Allreduce(&xi_ave, &xiave, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&xinormlocal, &xinorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  xinorm = sqrt(xinorm);
  
  FFTW_Slab_Conduct::update_conduct_order_param();
}

void FFTW_Slab_Conduct::update_vs()
{
  int lN1 = local_x;
  int N1 = nx;
  int N2 = local_y;
  int N3 = local_z;
  int NS = slip_systems;
  int ND = dimension;
  int index=0, na=0, na0=0, na1=0;
  double xinormlocal=0.0;           // stores the norm of the increment of xi
  double xirep=0.0, xiimp=0.0;      // store previous re and im values of xi
  double xi_ave=0.0;                // Local average of the order parameter
  int nsize = N1*N2*N3;
  double vs_diff_real, vs_diff_imag;
  int ida, idb, i, j, k, m, l, index2;
/* double mu= material->mu;

for(ida=0;ida<ND;ida++){
for (idb=0;idb<ND;idb++){
  for(i=0;i<lN1;i++){
    for(j=0;j<N2;j++){
      for (k=0;k<N3;k++){
        na0 = 2*(k + j*N3 + i*N2*N3 + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
        aveVs[ida][idb] += vstrain[na0]; //aveps only appears here to get total average strain in the box
      }
    }
  }
  MPI_Allreduce(&aveVs[ida][idb], &ave_Vs[ida][idb], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ave_Vs[ida][idb] /= ((double) nsize);
} //idb
} //ida
double sigma_Vs[3][3];
for (ida=0;ida<ND;ida++){
for (idb=0;idb<ND;idb++){
for (m=0;m<ND;m++){
for (l=0;l<ND;l++){
sigma_Vs[ida][idb] += C[ida][idb][m][l]*ave_Vs[m][l]/mu;
}
}
}
}*/

  for(ida=0;ida<ND;ida++) {
    for(idb=0;idb<ND;idb++){

      na = ida+idb*ND;

      for(i=0;i<lN1;i++) {
        for(j=0;j<N2;j++) {
          for(k=0;k<N3;k++) {

            na0 = 2*(i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
            index = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
            na1 = na0+1;
            index2 = i*N2*N3 + j*N3 + k;

            // elastic_energy_den_vs is the change made within one time step dtt

            if(phase[index2] != 0){

              vs_diff_real = (app->CD*app->timestep)*(elastic_energy_den_vs[index].re-sigma_rot[ida][idb]);//-sigma_Vs[ida][idb]);

              // vs_diff[na] += pow(vs_diff_real, 2.);

              vstrain[na0] = vstrain[na0]-vs_diff_real;

              vs_diff_imag = (app->CD*app->timestep)*elastic_energy_den_vs[index].im;

              vstrain[na1] = vstrain[na1]-vs_diff_imag;

            }

            // vstrain_complex[index].re = vstrain[na0];
            // vstrain_complex[index].im = vstrain[na1];

          }
        }
      } /*end i,j,k*/
    } // idb
  } // ida

}
/* -----------------------------------------------------------------------
   Update conductivity order parameter
   ---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::update_conduct_order_param()
{
  int lN1 = local_x;
  int N1 = nx;
  int N2 = local_y;
  int N3 = local_z;
  int index=0, na=0, na0=0, na1=0;
  double condxinormlocal=0.0;           // stores the norm of the increment of condxi
  double condxirep=0.0, condxiimp=0.0;      // store previous re and im values of condxi
  double condxi_ave=0.0;                // Local average of the conduct. order parameter
  int nsize = N1*N2*N3;

  //~ for(int isa=0;isa<NS;isa++){
	for(int i=0;i<lN1;i++){
	  for(int j=0;j<N2;j++){
	for(int k=0;k<N3;k++){
	  na0 = 2*(i*N2*N3 + j*N3 + k); // + isa*lN1*N2*N3);
	  index = i*N2*N3 + j*N3 + k; // + isa*lN1*N2*N3;
	  na1 = na0+1;

	  //Ginzburg-Landau Equation for real and imag parts

	  //if(xo[index] == 0.0){ // do we need this??
		condxirep = condxi[0][na0];
		condxiimp = condxi[0][na1];
		// assuming for now that -partial_i(sigma_ij E^ex_j)=0 and partial_i sigma_ij=0
		// and also take same time step size as for dislocation sub-problem (may generalize as needed)
		condxi[0][na0] = condxi[0][na0] -(1e-22*(app->timestep)*(conddata_real[index].re/(nsize))); // timestep is set to 1 in app, need smaller here
		condxi[0][na1] = condxi[0][na1] -(1e-22*(app->timestep)*(conddata_real[index].im/(nsize)));
	  //}
	  condxinormlocal += (condxi[0][na0] - condxirep)*(condxi[0][na0] - condxirep) +
		(condxi[0][na1] - condxiimp)*(condxi[0][na1] - condxiimp);
	  condxi_ave += condxi[0][na0];
	}
	  }
	}
  //~ }
  MPI_Allreduce(&condxi_ave, &condxiave, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&condxinormlocal, &condxinorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  condxinorm = sqrt(condxinorm);
}
/* ------------------------------------------------------------------
Prepares the next iteration
------------------------------------------------------------------ */
void FFTW_Slab_Conduct::prepare_next_itr()
{
  int lN1 = local_x;
  int N1 = nx;
  int N2 = local_y;
  int N3 = local_z;
  int NS = slip_systems;
  int ND = dimension;
  int index=0, na=0, na0=0, na1=0;
  int ida, idb, i, j, k;

  for(int isa=0;isa<NS;isa++){
    for(int i=0;i<lN1;i++){
      for(int j=0;j<N2;j++){
        for(int k=0;k<N3;k++){
          na0 = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          na1 = na0+1;
          data_fftw[index].re = xi[0][na0];
          data_fftw[index].im = xi[0][na1];

        }
      }
    }
  }
  if (app->nphase>1) {
    for(ida=0;ida<ND;ida++) {
      for(idb=0;idb<ND;idb++){
        for(i=0;i<lN1;i++) {
          for(j=0;j<N2;j++) {
            for(k=0;k<N3;k++) {
              na0 = 2*(i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
              index = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
              na1 = na0+1;
              vstrain_complex[index].re = vstrain[na0];
              vstrain_complex[index].im = vstrain[na1];

            }
          }
        }


      }
    }
  }
  FFTW_Slab_Conduct::prepare_next_conditr();
}
/* ------------------------------------------------------------------
   Calculate location of trailing partial for Peierls stress simulation
   ------------------------------------------------------------------ */
   void FFTW_Slab_Conduct::peierls(int p)
   {
     int lN1 = local_x;
     int lxs = local_x_start;
     int N1 = nx;
     int N2 = local_y;
     int N3 = local_z;
     int NS = slip_systems;
     int NP = num_planes;
     int ND = dimension;
     int ptindex=0, index=0, na=0, na0=0, na1=0, index1=0, index2=0, index3=0;
     int ida, idb, i, j, k, ip, isa, is;
     double peierlsdist = sqrt(3.0);

     int rank, np;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &np);
     MPI_Status status;
  
  if(solve->sstate == 0){
    if(p == 0){
      if(me == 0) printf("Looking for trailing partial...\n");
        for(ip=0;ip<NP;ip++){
          for(i=0;i<lN1;i++){
              for(j=0;j<N2;j++){
                  for(k=0;k<N3;k++){
                    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
                    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
                    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
                    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
                    na1=na0+1;

                    if(((lxs+i) == 0 && j == int(N2/2.0) && k < int(N3/2.0) && ip == 0) && fabs((sqrt(3.0)/2.0)-dregp[0][na0]) < peierlsdist){
                      peierlsz = k;
                      peierlsdreg=dregp[0][na0];
                      peierlsdist = fabs((sqrt(3.0)/2.0)-dregp[0][na0]);
                    }
                  }}}}
        if(me == 0) printf("Trailing partial found.\n");
        if(me == 0) printf("peierls z = %d\npeierls dreg = %f\n",peierlsz,peierlsdreg);
  }
    if(p == 1){
      if(me == 0) printf("Looking for leading partial...\n");
        for(ip=0;ip<NP;ip++){
          for(i=0;i<lN1;i++){
              for(j=0;j<N2;j++){
                  for(k=0;k<N3;k++){
                    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
                    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
                    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
                    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
                    na1=na0+1;

                    if(((lxs+i) == 0 && j == int(N2/2.0) && k < int(N3/2.0) && ip == 0) && fabs((sqrt(3.0)/2.0)-dregp[0][na1]) < peierlsdist){
                      peierlsz = k;
                      peierlsdreg=dregp[0][na1];
                      peierlsdist = fabs((sqrt(3.0)/2.0)-dregp[0][na1]);
                    }
                  }}}}
        if(me == 0) printf("Trailing partial found.\n");
        if(me == 0) printf("peierls z = %d\npeierls dreg = %f\n",peierlsz,peierlsdreg);
  }
  }
  }     

/* ------------------------------------------------------------------
   Calculate FCC screw dislocation core location
   ------------------------------------------------------------------ */
   void FFTW_Slab_Conduct::screw_center()
   {
     int lN1 = local_x;
     int lxs = local_x_start;
     int N1 = nx;
     int N2 = local_y;
     int N3 = local_z;
     int NS = slip_systems;
     int NP = num_planes;
     int ND = dimension;
     int ptindex=0, index=0, na=0, na0=0, na1=0, index1=0, index2=0, index3=0;
     int ida, idb, i, j, k, ip, isa, is;
     double peierlsdist = 2.0, centerdist1_1 = 1.0, centerdist1_2 = 1.0, centerdist2_1 = 1.0, centerdist2_2 = 1.0;
     int centerz1_mpi = 1, centerz2_mpi = 1;
     int centerx1_mpi = 0, centerx2_mpi = 0;
     int xslipflag = 0;

     int rank, np;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &np);
     MPI_Status status;

  // if(me == 0) printf("Looking for dislocation centers...\n");
  if(NP > 1){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
          for(k=0;k<N3;k++){
            index1=2*(i*N2*N3+j*N3+k+0*lN1*N2*N3);
            index2=2*(i*N2*N3+j*N3+k+1*lN1*N2*N3);
            index3=2*(i*N2*N3+j*N3+k+2*lN1*N2*N3);
            na0 = 2*(i*N2*N3+j*N3+k);
            na1=na0+1;

            if(fabs(dreg[0][index2]) >= 1){
              xslipflag = 1;
            }
          }}}
  }
  if(xslipflag == 1){
    for(i=0;i<lN1;i++){
      for(j=0;j<N2;j++){
          for(k=0;k<N3;k++){
            index1=2*(i*N2*N3+j*N3+k+0*lN1*N2*N3);
            index2=2*(i*N2*N3+j*N3+k+1*lN1*N2*N3);
            index3=2*(i*N2*N3+j*N3+k+2*lN1*N2*N3);
            na0 = 2*(i*N2*N3+j*N3+k);
            na1=na0+1;

            if(((lxs+i) == 0 && j == int(N2/2.0) && k < int(N3/2.0) && ip == 0) && fabs(1.0-dreg[0][na0]) < centerdist1_1){
              centerz1_mpi = k;
              centerdist1_1 = fabs(1.0-dreg[0][na0]);
            }
            if(((lxs+i) == 0 && j == int(N2/2.0) && k > int(N3/2.0) && ip == 0) && fabs(1.0-dreg[0][na0]) < centerdist1_2){
              centerz2_mpi = k;
              centerdist1_2 = fabs(1.0-dreg[0][na0]);
            }

            if(((lxs+i) > int(N1/2.0) && j == int(N2/2.0) && k == int(N3/8.0)) && fabs(1.0-fabs(dreg[0][index2])) < centerdist2_1){
              centerx1_mpi = lxs+i;
              centerdist2_1 = fabs(1.0-fabs(dreg[0][index2]));
            }
            if(((lxs+i) < int(N1/2.0) && j == int(N2/2.0) && k == int(N3*(7.0/8.0))) && fabs(1.0-fabs(dreg[0][index2])) < centerdist2_2){
              centerx2_mpi = lxs+i;
              centerdist2_2 = fabs(1.0-fabs(dreg[0][index2]));
            }
          }}}
  }
  else{
    for(i=0;i<lN1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
              index1=2*(i*N2*N3+j*N3+k+0*lN1*N2*N3);
              index2=2*(i*N2*N3+j*N3+k+1*lN1*N2*N3);
              index3=2*(i*N2*N3+j*N3+k+2*lN1*N2*N3);
              na0 = 2*(i*N2*N3+j*N3+k);
              na1=na0+1;

              if(((lxs+i) == 0 && j == int(N2/2.0) && k < int(N3/2.0)) && fabs(1.0-dreg[0][na0]) < centerdist1_1){
                centerz1_mpi = k;
                centerdist1_1 = fabs(1.0-dreg[0][na0]);
              }
              if(((lxs+i) == 0 && j == int(N2/2.0) && k > int(N3/2.0)) && fabs(1.0-dreg[0][na0]) < centerdist1_2){
                centerz2_mpi = k;
                centerdist1_2 = fabs(1.0-dreg[0][na0]);
              }
            }}}
          }
  // if(me == 0) printf("Dislocation centers found.\n");
  // if(me == 0) printf("Center z = %d,%d\n",centerz1_mpi,centerz2_mpi);

    if(np > 1){

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&centerz1, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&centerz1, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&centerz1_mpi, 1, MPI_INT, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&centerz1, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&centerz1, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&centerz1, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&centerz2, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&centerz2, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&centerz2_mpi, 1, MPI_INT, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&centerz2, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&centerz2, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&centerz2, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
   }
   if(np > 1){

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&centerx1, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&centerx1_mpi, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&centerx1, 1, MPI_INT, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&centerx1, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&centerx1, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&centerx1, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&centerx2, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&centerx2, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&centerx2_mpi, 1, MPI_INT, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&centerx2, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&centerx2, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&centerx2, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
   }
  }   

/* ------------------------------------------------------------------
   Initialize cross slip plane order parameters at dislocation centers (obstacle-free simulation)
   ------------------------------------------------------------------ */
void FFTW_Slab_Conduct::xslip_sxtal_6SFCCPlanar()
{
  int is, ip, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
  int na0, na, na1, ptindex, index, index1, index2, index3, indexm, *nodes, rtn_val0, rtn_val1, layer;
  double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
  FILE *of0;
  char infile[100], input[100], c[10];
  double xf=0, yf=0,zf=0;
  double eps=1e-1;

  int NS = slip_systems;
  int NP = num_planes;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  obsden = 0.1;

  if(me==0){
    if (logfile) fprintf(logfile,"Initializing cross slip planes ...\n");
    if (screen) fprintf(screen,"Initializing cross slip planes ...\n");
  }
  for(is=0;is<NS;is++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na = 2*(i*N2*N3+j*N3+k);
    index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
    ptindex = i*N2*N3+j*N3+k;
    na1 = na0+1;

    if(is < 3 && (lxs+i) == 0 /*&& ((k > int(N3/8.0)) && k < int(N3*7.0/8.0))*/){
      xo[index] = 0.0;
    }
    else if (is > 2 && ((k == centerz1 && k != 0) || (k == centerz2 && k != 0))){
      xo[index] = 0.0;
    }
    else{
      xo[index] = 1.0;
    }
    // printf("opc[%d][%d] VS %d\n",is,ptindex,index);
    // opc[is][ptindex] = std::complex<double> (xi[0][na0],xi[0][na1]);
  //   data_fftw[index].re = xi[0][na0];
  //   data_fftw[index].im = xi[0][na1];
  }
  if(me==0){
    if (logfile) fprintf(logfile,"Initial configuration of cross slip planes COMPLETE ...\n");
    if (screen) fprintf(screen,"Initial configuration of cross slip planes COMPLETE  ...\n");
  }
  return;
} 

/* ------------------------------------------------------------------
   Check for screw annihilation/cross slip
   ------------------------------------------------------------------ */
   void FFTW_Slab_Conduct::xslip_check()
   {
     int lN1 = local_x;
     int lxs = local_x_start;
     int N1 = nx;
     int N2 = local_y;
     int N3 = local_z;
     int NS = slip_systems;
     int NP = num_planes;
     int ND = dimension;
     int ptindex=0, index=0, na=0, na0=0, na1=0, index1=0, index2=0, index3=0;
     int ida, idb, i, j, k, ip, isa, is;
     double Hflag_mpi = 0.0, X1flag_mpi = 0.0, X2flag_mpi = 0.0;
    //  Hflag = 2;

     int rank, np;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &np);
     MPI_Status status;

  for(ip=0;ip<NP;ip++){
    for(i=0;i<lN1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
              index = i*N2*N3+j*N3+k;
              index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
              index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
              index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
              na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
              na1=na0+1;

              if(((lxs+i) == 0 && (k == int(N3/8.0)-1 || k == int(N3*7.0/8.0)+1) && ip == 0)){
              // if(((lxs+i) == 0 && (k == 0 || k == N3-1) && ip == 0)){
              // if(((lxs+i) == 0 && (k == int(N3/4.0)-10 || k == int(N3*3.0/4.0)+10) && ip == 0)){ // Peierls stress/Timestep testing
              // if(((lxs+i) == 0 && (k == int(N3*5.0/16.0)-1 || k == int(N3*11.0/16.0)+1) && ip == 0)){ // 64x64x64 --> 64x64x128
                if(abs(dreg[0][na0]) >= 2.0) Hflag_mpi += 2.0;
                else Hflag_mpi += abs(dreg[0][na0]);
              }
              // if(((lxs+i) == 0 && (k == int(N3*5.0/16.0)-1 || k == int(N3*11.0/16.0)+1) && ip == 0)) Hflag += abs(dreg[0][na0]); // 64x64x128 --> 128x64x256
              // if(((lxs+i) == 0 && (k == int(N3*9.0/40.0) || k == int(N3*31.0/40.0)) && ip == 0)) Hflag += abs(dreg[0][na0]);
              // if(((lxs+i) == 0 && (k == N3-1 || k == 0) && ip == 0)) Hflag += abs(dreg[0][na0]);
              // if(((lxs+i) == 0 && (k == centerz1-1 || k == centerz2+1) && ip == 0)) Hflag += abs(dreg[0][na0]);
              if((((lxs+i) == N1-1/*int(N1/np)*/ && (k == int(N3/8.0) && k != 0))) && ip == 1){
              // if((((lxs+i) == 1/*int(N1/np)*/ && (k == int(N3/8.0) && k != 0))) && ip == 1){
              // if((((lxs+i) == N1-1/*int(N1/np)*/ && (k == int(N3*5.0/16.0) && k != 0))) && ip == 1){ // 64x64x64 --> 64x64x128
              // if((((lxs+i) == int(N1/np)-1 && k == int(N3*7.0/8.0))) && ip == 1){
              if(abs(dreg[0][na0]) >= 2.0) X1flag_mpi += 2.0;
              else X1flag_mpi += abs(dreg[0][na0]);
              // if(abs(dreg[0][na0]) < 2) X1flag_mpi = 0;
              }
              if((((lxs+i) == 1/*int(N1/np)-1*/ && (k == int(N3*7.0/8.0) && k != 0))) && ip == 1){
              // if((((lxs+i) == N1-1/*int(N1/np)-1*/ && (k == int(N3*7.0/8.0) && k != 0))) && ip == 1){
              // if((((lxs+i) == 1/*int(N1/np)-1*/ && (k == int(N3*(11.0/16.0)) && k != 0))) && ip == 1){ // 64x64x64 --> 64x64x128
              // if((((lxs+i) == N1-int(N1/np) && k == int(N3/8.0))) && ip == 1){
              if(abs(dreg[0][na0]) >= 2.0) X2flag_mpi += 2.0;
              else X2flag_mpi += abs(dreg[0][na0]);
              // if(abs(dreg[0][na0]) < 2) X2flag_mpi = 0;
              }
            }}}}
    if(np > 1){

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&X1flag, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&X1flag_mpi, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&X1flag, 1, MPI_DOUBLE, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&X1flag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    // else if(rank == int(np/4)-1){
    //   MPI_Recv(&X1flag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
    //   MPI_Send(&X1flag_mpi, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    // }
    else{
      MPI_Recv(&X1flag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&X1flag, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&X2flag, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&X2flag, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&X2flag_mpi, 1, MPI_DOUBLE, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&X2flag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    // else if(rank == int(np*3/4)){
    //   MPI_Recv(&X2flag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
    //   MPI_Send(&X2flag, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    // }
    else{
      MPI_Recv(&X2flag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&X2flag, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&Hflag, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&Hflag, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&Hflag_mpi, 1, MPI_DOUBLE, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&Hflag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    // else if(rank == int(np*3/4)){
    //   MPI_Recv(&X2flag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
    //   MPI_Send(&X2flag, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    // }
    else{
      MPI_Recv(&Hflag, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&Hflag, 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
   }
    Hflag = Hflag/(2*N2);
    X1flag = X1flag/(N2);
    X2flag = X2flag/(N2);

  }

/* ------------------------------------------------------------------
   Copy barrier xo to all slabs for kMC simulations
   ------------------------------------------------------------------ */
   void FFTW_Slab_Conduct::barrier_xslip()
   {
     int lN1 = local_x;
     int lxs = local_x_start;
     int N1 = nx;
     int N2 = local_y;
     int N3 = local_z;
     int NS = slip_systems;
     int NP = num_planes;
     int ND = dimension;
     int ptindex=0, index=0, na=0, na0=0, na1=0, index1=0, index2=0, index3=0, index4=0, index5=0, index6=0;
     int ida, idb, i, j, k, ip, isa, is;
     int nmpi = N2*N3;
     int bar[nmpi], bar_mpi[nmpi], bar2[nmpi], bar2_mpi[nmpi];

     for(i=0;i<nmpi;i++){
      bar[i] = 0.0;
      bar_mpi[i] = 0.0;
     }
     is = 0;
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
    index1 = i*N2*N3 + j*N3 + k + (is+0)*lN1*N2*N3;
    index2 = i*N2*N3 + j*N3 + k + (is+1)*lN1*N2*N3;
    index3 = i*N2*N3 + j*N3 + k + (is+2)*lN1*N2*N3;
    index4 = i*N2*N3 + j*N3 + k + (is+3)*lN1*N2*N3;
    index5 = i*N2*N3 + j*N3 + k + (is+4)*lN1*N2*N3;
    index6 = i*N2*N3 + j*N3 + k + (is+5)*lN1*N2*N3;
    index = j*N3 + k;

    if(xo[index4] == 1.0 && dt_kmc[index2] > 1.0) bar_mpi[index] = 1.0;
    if(xo[index4] == -1.0) bar_mpi[index] = -1.0;
    if(xo[index4] == 2.0){
      bar_mpi[index] = 2.0;
      bar2_mpi[index] = 2.0;
    }
    }

     int rank, np;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &np);
     MPI_Status status;

    if (rank == np-1) {
      MPI_Recv(&bar, nmpi, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&bar, nmpi, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&bar_mpi, nmpi, MPI_INT, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&bar, nmpi, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&bar, nmpi, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&bar, nmpi, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == np-1) {
      MPI_Recv(&bar2, nmpi, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&bar2_mpi, nmpi, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    else if(rank == 0){
      MPI_Send(&bar2, nmpi, MPI_INT, np-1, 1, MPI_COMM_WORLD);
      MPI_Recv(&bar2, nmpi, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
    }
    else{
      MPI_Recv(&bar2, nmpi, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
      MPI_Send(&bar2, nmpi, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1 = i*N2*N3 + j*N3 + k + (is+0)*lN1*N2*N3;
    index2 = i*N2*N3 + j*N3 + k + (is+1)*lN1*N2*N3;
    index3 = i*N2*N3 + j*N3 + k + (is+2)*lN1*N2*N3;
    index4 = i*N2*N3 + j*N3 + k + (is+3)*lN1*N2*N3;
    index5 = i*N2*N3 + j*N3 + k + (is+4)*lN1*N2*N3;
    index6 = i*N2*N3 + j*N3 + k + (is+5)*lN1*N2*N3;
    index = j*N3 + k;

    if(bar[index] == 1.0){
      xo[index4] = 1.0;
      xo[index5] = 1.0;
      xo[index6] = 1.0;
    }
    if(bar[index] == -1.0){
      xo[index4] = -1.0;
      xo[index5] = -1.0;
      xo[index6] = -1.0;
    }
    if(bar[index] == 2.0 || bar2[index] == 2.0){
      xo[index4] = 2.0;
      xo[index5] = 2.0;
      xo[index6] = 2.0;
    }
  }
  }

		/* ------------------------------------------------------------------
		   Prepares the next conductivity iteration
		   ------------------------------------------------------------------ */
		void FFTW_Slab_Conduct::prepare_next_conditr()
		{
		  int lN1 = local_x;
		  int N1 = nx;
		  int N2 = local_y;
		  int N3 = local_z;
		  int index=0, na=0, na0=0, na1=0;

			for(int i=0;i<lN1;i++){
			  for(int j=0;j<N2;j++){
			for(int k=0;k<N3;k++){
			  na0 = 2*(i*N2*N3 + j*N3 + k);
			  index = i*N2*N3 + j*N3 + k;
			  //~ na = 2*(i*N2*N3 + j*N3 + k);
			  na1 = na0+1;
			  conddata_fftw[index].re = condxi[0][na0];
			  conddata_fftw[index].im = condxi[0][na1];

			  }
			}
		  }
		}
        /* -----------------------------------------------------------------------
        temporary stores data_fftw
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::temp_stor_data()
        {
          int NS = slip_systems;
          int index = -1;
          for(int isa=0;isa<NS;isa++){
            for(int ii=0;ii<total_local_size;ii++){
              index = ii + isa*total_local_size;
              temp_data[index].re = data_fftw[index].re;
              temp_data[index].im = data_fftw[index].im;
            }
          }
          if (fft->me == 0) {
            if (screen) fprintf(screen,"data_fftw->temp_data ...\n");
            if (logfile) fprintf(logfile,"data_fftw->temp_data ...\n");
          }
        }
        /* -----------------------------------------------------------------------
        reset data_fftw
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::reset_data()
        {
          int NS = slip_systems;
          int index = -1;
          for(int isa=0;isa<NS;isa++){
            for(int ii=0;ii<total_local_size;ii++){
              index = ii + isa*total_local_size;
              data_fftw[index].re = temp_data[index].re;
              data_fftw[index].im = temp_data[index].im;
            }
          }
          if (fft->me == 0) {
            if (screen) fprintf(screen,"temp_data->data_fftw ...\n");
            if (logfile) fprintf(logfile,"temp_data->data_fftw ...\n");
          }
        }

        /* ------------------------------------------------------------------
        Prepares the FFT forward
        ------------------------------------------------------------------ */
        void FFTW_Slab_Conduct::prep_forward()
        {
          if(mode == 1)
          forward_mode1();
          else if(mode == 2)
          forward_mode2();
          FFTW_Slab_Conduct::forwardcond();
        }

        /* ------------------------------------------------------------------
        Prepares the FFT backward
        ------------------------------------------------------------------ */
        void FFTW_Slab_Conduct::prep_backward()
        {
          if(mode == 1)
          backward_mode1();
          else if(mode == 2)
          backward_mode2();
          FFTW_Slab_Conduct::backwardcond();
        }


        /* -----------------------------------------------------------------------
        forward FFT
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::forward_mode1()
        {
          int lN1 = local_x;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;

          for(int i=0; i<NS; i++){
            int psys = i*lN1*N2*N3;
            fftwnd_mpi(plan, 1, data_fftw+psys, work, FFTW_NORMAL_ORDER);
          }
        }

        /* -----------------------------------------------------------------------
        forward FFT
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::forward_mode2()
        {
          int lN1 = local_x;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;
          int index = -1, index2 = -1;


          for(int isa=0;isa<NS;isa++){
            for(int ii=0;ii<total_local_size;ii++){
              index = isa + ii*NS;
              index2 = ii + isa*total_local_size;
              temp_data[index] = data_fftw[index2];
            }
          }
          fftwnd_mpi(plan,NS,temp_data,work,FFTW_NORMAL_ORDER);

          for(int isa=0;isa<NS;isa++){
            for(int ii=0;ii<total_local_size;ii++){
              index = isa + ii*NS;
              index2 = ii + isa*total_local_size;
              data_fftw[index2] = temp_data[index];
            }
          }
        }

		/* -----------------------------------------------------------------------
		   forward FFT - conductivity order parameter
		   ---------------------------------------------------------------------*/
		void FFTW_Slab_Conduct::forwardcond()
		{
		  int lN1 = local_x;
		  int N2 = local_y;
		  int N3 = local_z;

		  fftwnd_mpi(plan, 1, conddata_fftw, work, FFTW_NORMAL_ORDER);
		  //~ fftwnd_mpi(plan, 1, jcurrent, work, FFTW_NORMAL_ORDER); // don't need this
		}

        /* -----------------------------------------------------------------------
        backward FFT
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::backward_mode1()
        {
          int lN1 = local_x;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;

          for(int i=0; i<NS; i++){
            int psys = i*lN1*N2*N3;
            fftwnd_mpi(iplan, 1, data_real+psys, work, FFTW_NORMAL_ORDER);
          }
        }

        /* -----------------------------------------------------------------------
        backward FFT
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::backward_mode2()
        {
          int lN1 = local_x;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;
          int index = -1, index2 = -1;

          for(int isa=0;isa<NS;isa++){
            for(int ii=0;ii<total_local_size;ii++){
              index = isa + ii*NS;
              index2 = ii + isa*total_local_size;
              temp_data[index] = data_real[index2];
            }
          }
          fftwnd_mpi(iplan,NS,temp_data,work,FFTW_NORMAL_ORDER);

          for(int isa=0;isa<NS;isa++){
            for(int ii=0;ii<total_local_size;ii++){
              index = isa + ii*NS;
              index2 = ii + isa*total_local_size;
              data_real[index2] = temp_data[index];
            }
          }
        }

		/* -----------------------------------------------------------------------
		   backward FFT - conductivity order parameter, and assemble parts in real space
		   ---------------------------------------------------------------------*/
		void FFTW_Slab_Conduct::backwardcond()
		{
		  int lN1 = local_x;
		  int N2 = local_y;
		  int N3 = local_z;
		  int index;
		  int nsize = lN1*N2*N3;
		  double epsilonzero = material->epsilonzero;
		  double condsigma=0, rhointerface=0, rhograin=0, rhodislocation=0;
		  double condsigma_array[N2];
		  double lam = 0.; // mean free path of electrons
		  //~ double d0 = 50e-9; // length scale, e.g. average layer thickness
		  double Vf = 0.; // volume fraction of materialB (calculated below)
		  int Nlayers;
		  int NgrainsA, NgrainsB, Ngrains;
		  int dlayerA, dlayerB;
		  double prob = 0.5; // one minus scattering probability, (1-p) (overwritten by user input)
		  double PSprob = 0;
		  double alphagrain = 0.;
		  double dT = app->temperature-300;
		  //
		  int rank, np, slabsize;
		  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		  MPI_Comm_size(MPI_COMM_WORLD, &np);
		  //~ slabsize = np*lN1;
		  
		  if (Nmaterials < 2) error->all(FLERR,"Need two materials for this app!");
		  // always have an equal number of A / B layers that alternate starting with material A:
		  Vf = (materialB->dlayer) / ((materialA->dlayer) + (materialB->dlayer));
		  Nlayers = materialA->Nlayers;
		  prob = 1.0-materialA->pinterface;
		  NgrainsA = materialA->Ngrains;
		  NgrainsB = materialB->Ngrains;
		  // layer thicknesses in integer (corresponding to number of x-indices per layer), rounded: 
		  dlayerB = int(round((N2/(Nlayers))*Vf));
		  dlayerA = int(round((N2/(Nlayers))*(1.-Vf)));
		  if (rank==0) printf("layermath: %i, %i\n",dlayerA,dlayerB);
		  if (rank==0) printf("Ngrains: %i, %i\n",NgrainsA,NgrainsB);

		  fftwnd_mpi(iplan, 1, conddata_bpart, work, FFTW_NORMAL_ORDER); // b_{ij} in notes
		  fftwnd_mpi(iplan, 1, jcurrent_fftw, work, FFTW_NORMAL_ORDER); // a_j in notes
		  for (int i=0;i<dimension;i++)
		  {
			javerage[i].re = 0;
			javerage[i].im = 0;
			jcurrent[i].re = 0;
			jcurrent[i].im = 0;
			conddata_real[i].re = 0.;
			conddata_real[i].im = 0.;
		  }
		  if (rank==0) printf("\ntemperature: %f \n", dT+300);
		  if (rank==0) printf("\nsize in y: %i Nmaterials: %i \n", N2, Nmaterials);
		  if (rank==0) printf("materialA: conductivity = %f MS/m, ", materialA->condsigma/1e6);
		  if (rank==0) printf("mean-free-path = %f nm, dlayer = %f nm\n", 1.e9*materialA->mfpel, 1.e9*materialA->dlayer);
		  if (rank==0) printf("materialB: conductivity = %f MS/m, ", materialB->condsigma/1e6);
		  if (rank==0) printf("mean-free-path = %f nm, dlayer = %f nm\n", 1.e9*materialB->mfpel, 1.e9*materialB->dlayer);
		  if (rank==0) printf("Nlayers: %i, Vf = %f \n",Nlayers,Vf);
		  printf("# cpus: %i, current process #: %i\n",np,rank);
		  
		  for(int j=0;j<N2;j++){
			if ((j)  % (dlayerA+dlayerB) >= dlayerA)
			{
			  //~ printf("Nb: %i, %i, %i, %i, %i",rank,j,(rank*np+i)  % (dlayerA+dlayerB), dlayerA, dlayerB);
			  condsigma_array[j] = materialB->condsigma / (1. + materialB->alpha_resistivity*dT);
			  PSprob = 2*prob/(materialB->dlayer);
			  lam = (materialB->mfpel)/(1. + materialB->alpha_mfpel*dT); // mean free path for Nb with temperature correction
			  // and with correction for misaligned layers (vary angles < 1/2 ~ 29 degrees):
			  //lam = (materialB->mfpel)*(1. + exp(-(materialB->dlayer)/abs(cos(0.5*j/lN1))/(materialB->mfpel)))/(1. + materialB->alpha_mfpel*dT);
			  alphagrain=(lam*materialB->Ngrains/materialB->dlayer)*materialB->Rgrain/(1-materialB->Rgrain);
			  Ngrains = materialB->Ngrains;
			  rhodislocation = materialB->beta_disloc*materialB->rho_disloc;
			  
			}
			else 
			{
			  condsigma_array[j] = materialA->condsigma / (1. + materialA->alpha_resistivity*dT);
			  //PSprob = 4*prob*Vf/((materialA->dlayer)*(1-Vf));
			  PSprob = 2*prob/(materialA->dlayer);
			  lam = (materialA->mfpel)/(1. + materialA->alpha_mfpel*dT); // mean free path of electrons in copper with temperature correction
			  // and with correction for misaligned layers (vary angles < 1/2 ~ 29 degrees):
			  //lam = (materialA->mfpel)*(1. + exp(-(materialA->dlayer)/abs(cos(0.5*j/lN1))/(materialA->mfpel)))/(1. + materialA->alpha_mfpel*dT);
			  alphagrain=(lam*materialA->Ngrains/materialA->dlayer)*materialA->Rgrain/(1-materialA->Rgrain);
			  Ngrains = materialA->Ngrains;
			  rhodislocation = materialA->beta_disloc*materialA->rho_disloc;
			}
			if (PSprob > 0.)
			{
			  //condsigma = 1/((1/condsigma_array[j]) *(1 + (3./16.)*lam*PSprob));
			  rhointerface = ((1/condsigma_array[j]) *( (3./16.)*lam*PSprob));
			  if (Ngrains > 1)
			  {
				rhograin = (1/condsigma_array[j]) *(1./(1.-3.*alphagrain/2.+3.*alphagrain*alphagrain - 3.*pow(alphagrain,3.)*log(1.+1./alphagrain)) - 1.);
				// weight the two effects by the number of interfaces occuring in the region; also add resistivity from dislocation density
				condsigma = 1/((1/condsigma_array[j]) + (2*rhointerface + (Ngrains-1)*rhograin)/(Ngrains+1) + rhodislocation);
			  }
			  else
			  {
				rhograin=0;
				condsigma = 1/((1/condsigma_array[j]) + rhointerface + rhodislocation);
			  }
			//		printf("\nlocalcond at x= %i: %f",j,condsigma);
			}
			else condsigma = condsigma_array[j]; // TODO: make this a tensor
			for(int i=0;i<lN1;i++){
			  for(int k=0;k<N3;k++){
				index  = i*N2*N3 + j*N3 + k;
				//~ for (int l=0;l<dimension;l++){
				jcurrent[index].re = condsigma*(material->El[1] - jcurrent_fftw[index].re/epsilonzero); // TODO: make this a vector
				jcurrent[index].im = - condsigma*jcurrent_fftw[index].im/epsilonzero;
				//~ }
				javerage[1].re += jcurrent[index].re/nsize; // average current needed for macroscopic conductivity
				javerage[1].im += jcurrent[index].im/nsize;
				conddata_real[index].re = (condsigma/epsilonzero)*conddata_bpart[index].re;
				conddata_real[index].im = (condsigma/epsilonzero)*conddata_bpart[index].im;
				if (j==0) {
					conddata_real[index].re += (jcurrent_fftw[index].re/epsilonzero - material->El[1])*(condsigma_array[j+1]-condsigma_array[j]); // TODO: extrapolate properly at the edge
					conddata_real[index].im += (jcurrent_fftw[index].im/epsilonzero)*(condsigma_array[j+1]-condsigma_array[j]);
				}
				else if (j==N2-1) {
					conddata_real[index].re += (jcurrent_fftw[index].re/epsilonzero - material->El[1])*(condsigma_array[j]-condsigma_array[j-1]); // TODO: extrapolate properly at the edge
					conddata_real[index].im += (jcurrent_fftw[index].im/epsilonzero)*(condsigma_array[j]-condsigma_array[j-1]);
				}
				else {
					conddata_real[index].re += (jcurrent_fftw[index].re/epsilonzero - material->El[1])*(condsigma_array[j+1]-condsigma_array[j-1])/2.; // TODO: divide by delta x: find out which variable stores positions!
					conddata_real[index].im += (jcurrent_fftw[index].im/epsilonzero)*(condsigma_array[j+1]-condsigma_array[j-1])/2.;
				}
			}}}
		   
		   conductivity = 0.;
		   for (int l=1;l<2;l++){ // TODO: update for vector
			 conductivity +=javerage[l].re / (double(np)*material->El[l]); // determine macroscopic conductivity per cpu (TODO: make this a tensor)
		   }
		   //~ conductivity = conductivity/((double)dimension);
		   MPI_Allreduce(&conductivity, &conductivity_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		}

        /* -----------------------------------------------------------------------
        Initial Crystal to nothing
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::initial_sxtal_NoLoop()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
              na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
              na = 2*(i*N2*N3+j*N3+k);
              index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
              na1 = na0+1;
              xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
              xi[0][na1] = 0.0;
              xo[index] = 0.0;
              data_fftw[index].re = xi[0][na0];
              data_fftw[index].im = xi[0][na1];
            }
          }

          if (app->nphase == 2) set_phase_2Px();
          if (app->nphase == 3) set_phase_3Px();

          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }

        /* -----------------------------------------------------------------------
        Initial Crystal NonorthogonalFCC/BCC 2 Loops in k and j plane
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::initial_sxtal_mode1()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;
          double ar[3];

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
              for(j=0;j<N2;j++)
        	for(k=0;k<N3;k++){
        	  na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
        	  na = 2*(i*N2*N3+j*N3+k);
        	  index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
        	  na1 = na0+1;
        	  xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
        	  xi[0][na1] = 0.0;
        	  xo[index] = 0.0;
        	  data_fftw[index].re = 0.0;
        	  data_fftw[index].im = 0.0;


        	  if(((k-(local_z/2.0))<eps) &&     // {111} plane
        	     ((k-(local_z/2.0))>-eps) &&
        	     (is==0)){

              ar[0] = sclprim[0][0]*((lxs+i)-N1/2.0)+sclprim[0][1]*(j-5.0*local_y/8.0)+sclprim[0][2]*(k-local_z/2.0);
              ar[1] = sclprim[1][0]*((lxs+i)-N1/2.0)+sclprim[1][1]*(j-5.0*local_y/8.0)+sclprim[1][2]*(k-local_z/2.0);
              ar[2] = sclprim[2][0]*((lxs+i)-N1/2.0)+sclprim[2][1]*(j-5.0*local_y/8.0)+sclprim[2][2]*(k-local_z/2.0);

              ir = sqrt(ar[0]*ar[0]+ar[1]*ar[1]+ar[2]*ar[2]);
        	    // ir = sqrt(((lxs+i)-N1/2.0)*((lxs+i)-N1/2.0)+(j-3.0*local_y/4.0)*(j-3.0*local_y/4.0)+(k-local_z/2.0)*(k-local_z/2.0));

        	    if(ir<=(local_y/8.0-2.0)){
        	      xi[0][na0]=1.0;
        	      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
        	    }
        	    //write input file
        	    //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
        	    data_fftw[index].re = xi[0][na0];
        	    data_fftw[index].im = xi[0][na1];
        	  }
        	  else if(((j-local_y/2.0)<eps) &&                              // {100} plane
        		  ((j-local_y/2.0)>-eps) &&
        		  (is==1)){
                ar[0] = sclprim[0][0]*((lxs+i)-N1/2.0+N1/8.0)+sclprim[0][1]*(j-local_y/2.0)+sclprim[0][2]*(k-5.0*local_z/8.0);
                ar[1] = sclprim[1][0]*((lxs+i)-N1/2.0+N1/8.0)+sclprim[1][1]*(j-local_y/2.0)+sclprim[1][2]*(k-5.0*local_z/8.0);
                ar[2] = sclprim[2][0]*((lxs+i)-N1/2.0+N1/8.0)+sclprim[2][1]*(j-local_y/2.0)+sclprim[2][2]*(k-5.0*local_z/8.0);

                ir = sqrt(ar[0]*ar[0]+ar[1]*ar[1]+ar[2]*ar[2]);

        	    // ir = sqrt(((lxs+i)-N1/2.0)*((lxs+i)-N1/2.0)+(j-local_y/2.0)*(j-local_y/2.0)+(k-3.0(k-3.0*local_z/8.0));
        	    if(ir<=(local_z/8.0-2.0)){
        	      xi[0][na0]=1.0;
        	      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
        	    }
        	    //write input file
        	    //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
        	    data_fftw[index].re = xi[0][na0];
        	    data_fftw[index].im = xi[0][na1];
        	  }
        	}
          }
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration MODE1 of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration Mode1 of order parameters COMPLETE  ...\n");
          }
          return;
        }

        /* -----------------------------------------------------------------------
        Initial Crystal mode2
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::initial_sxtal_mode2()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;
          // double ar[3];

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
              for(j=0;j<N2;j++)
        	for(k=0;k<N3;k++){
        	  na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
        	  na = 2*(i*N2*N3+j*N3+k);
        	  index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
        	  na1 = na0+1;
        	  xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
        	  xi[0][na1] = 0.0;
        	  xo[index] = 0.0;
        	  data_fftw[index].re = 0.0;
        	  data_fftw[index].im = 0.0;


            if (((lxs+i+j+k-N1)<eps) && ((lxs+i+j+k-N1)>-eps) &&
            (is==0)){
              ir = app->sa*sqrt(((lxs+i)-28.0)*((lxs+i)-28.0)+(j-28.0)*(j-28.0)+(k-40.0)*(k-40.0));
              if(ir <= 3.0*sqrt(2.0)){
        	      xi[0][na0]=1.0;
        	      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
        	    }
        	    //write input file
        	    //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
        	    data_fftw[index].re = xi[0][na0];
        	    data_fftw[index].im = xi[0][na1];
        	  }
        	  else if (((lxs+i+j-k)<eps) && ((lxs+i+j-k)>-eps) &&
        	  (is==1)){
              ir = app->sa*sqrt(((lxs+i)-28.0)*((lxs+i)-28.0)+(j-28.0)*(j-28.0)+(k-56.0)*(k-56.0));
              if(ir <= 3.0*sqrt(2.0)){
        	      xi[0][na0]=1.0;
        	      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
        	    }
        	    //write input file
        	    //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
        	    data_fftw[index].re = xi[0][na0];
        	    data_fftw[index].im = xi[0][na1];
        	  }
        	}
          }
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration MODE1 of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration Mode1 of order parameters COMPLETE  ...\n");
          }
          return;
        }

        // This function generates 1 shear loop in an orthogonal grid
        //  oriented with [1-10], [11-2] and [111] axis for fcc
        //  so the loop gliding plane will be the {111}w

        void FFTW_Slab_Conduct::initial_sxtal_1LOrthoFCC()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
              na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
              na = 2*(i*N2*N3+j*N3+k);
              index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
              na1 = na0+1;
              xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
              xi[0][na1] = 0.0;
              xo[index] = 0.0;
              data_fftw[index].re = 0.0;
              data_fftw[index].im = 0.0;

              if(((k-(local_z/2.0))<eps) &&                              // {100} plane
              ((k-(local_z/2.0))>-eps) &&
              (is==0)){
                ir = sqrt(((lxs+i)-N1/2.0)*((lxs+i)-N1/2.0)+(j-local_y/2.0)*(j-local_y/2.0));

                if(ir<=xprd/10.0){
                  xi[0][na0]=1.0;
                  xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
                }
                //write input file
                //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
                data_fftw[index].re = xi[0][na0];
                data_fftw[index].im = xi[0][na1];
              }
            }
          }
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }
        // A loop on (1 1 1) in non-orthogonal FCC grid
        void FFTW_Slab_Conduct::initial_sxtal_1LNonOrthoFCC()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;
          double ar[3];

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
              for(j=0;j<N2;j++)
        	for(k=0;k<N3;k++){
        	  na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
        	  na = 2*(i*N2*N3+j*N3+k);
        	  index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
        	  na1 = na0+1;
        	  xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
        	  xi[0][na1] = 0.0;
        	  xo[index] = 0.0;
        	  data_fftw[index].re = 0.0;
        	  data_fftw[index].im = 0.0;

            if(((k-(local_z/2.0))<eps) &&
            ((k-(local_z/2.0))>-eps) &&
            (is==0)){
              ar[0] = sclprim[0][0]*((lxs+i)-N1/2.0)+sclprim[0][1]*(j-local_y/2.0)+sclprim[0][2]*(k-local_z/2.0);
              ar[1] = sclprim[1][0]*((lxs+i)-N1/2.0)+sclprim[1][1]*(j-local_y/2.0)+sclprim[1][2]*(k-local_z/2.0);
              ar[2] = sclprim[2][0]*((lxs+i)-N1/2.0)+sclprim[2][1]*(j-local_y/2.0)+sclprim[2][2]*(k-local_z/2.0);

              ir = sqrt(ar[0]*ar[0]+ar[1]*ar[1]+ar[2]*ar[2]);
              if(ir <= 10.0){
                xi[0][na0]=1.0;
                xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
              }
              data_fftw[index].re = xi[0][na0];
              data_fftw[index].im = xi[0][na1];
            }
        	} // i, j, k
        } // is
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }

        /* -----------------------------------------------------------------------
        Initial Crystal fcc with 1 shear loop on the inclined (1 1 1) plane
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::initial_sxtal_1LInclinedFCC()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
              for(j=0;j<N2;j++)
        	for(k=0;k<N3;k++){
        	  na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
        	  na = 2*(i*N2*N3+j*N3+k);
        	  index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
        	  na1 = na0+1;
        	  xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
        	  xi[0][na1] = 0.0;
        	  xo[index] = 0.0;
        	  data_fftw[index].re = 0.0;
        	  data_fftw[index].im = 0.0;

            if (((lxs+i+j+k-N1)<eps) && ((lxs+i+j+k-N1)>-eps) &&
            (is==0)){
              ir = app->sa*sqrt(((lxs+i)-N1/3.0)*((lxs+i)-N1/3.0)+(j-local_y/3.0)*(j-local_y/3.0)+(k-local_z/3.0)*(k-local_z/3.0));
              if(ir <= 5.0*sqrt(2.0)){
        	      xi[0][na0]=1.0;
        	      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
        	    }
        	    data_fftw[index].re = xi[0][na0];
        	    data_fftw[index].im = xi[0][na1];
            }
        	} //i,j,k
        } //NS
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }
        /* -----------------------------------------------------------------------
        Initial Crystal bcc with 1 shear loop in a crystal with orientation
        [1 1 -2], [1 1 1] and [1-10]
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::initial_sxtal_1LOrthoBCC()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
              na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
              na = 2*(i*N2*N3+j*N3+k);
              index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
              na1 = na0+1;
              xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
              xi[0][na1] = 0.0;
              xo[index] = 0.0;
              data_fftw[index].re = 0.0;
              data_fftw[index].im = 0.0;

              if(((k-(local_z/2.0))<eps) &&                              // {100} plane
              ((k-(local_z/2.0))>-eps) &&
              (is==0)){
                ir = sqrt(app->sa*((lxs+i)-N1/2.0)*app->sa*((lxs+i)-N1/2.0)+app->sb*(j-local_y/2.0)*app->sb*(j-local_y/2.0));

                //if(ir<=xprd/6.0){
                if(ir <= 10*app->sa){
                  xi[0][na0]=1.0;
                  xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
                }
                //write input file
                //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
                data_fftw[index].re = xi[0][na0];
                data_fftw[index].im = xi[0][na1];
              }
            }
          }
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }

        /* -----------------------------------------------------------------------
        Initial Crystal mode1
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::initial_sxtal_2LOrtho()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
              na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
              na = 2*(i*N2*N3+j*N3+k);
              index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
              na1 = na0+1;
              xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
              xi[0][na1] = 0.0;
              xo[index] = 0.0;

              // Rotate to grid
              // xf = ((lxs+i)*four[0][0] + j*four[0][1] + k*four[0][2])/normfour[0];
              // yf = ((lxs+i)*four[1][0] + j*four[1][1] + k*four[1][2])/normfour[1];
              // zf = ((lxs+i)*four[2][0] + j*four[2][1] + k*four[2][2])/normfour[2];

              if((2.8284*(j-local_y/2.0)+(k-local_z/2.0)<20.0*eps) &&     // {111} plane
              (2.8284*(j-local_y/2.0)+(k-local_z/2.0)>-20.0*eps) &&
              (is==1)){
                ir = sqrt(((lxs+i)-N1/2.0)*((lxs+i)-N1/2.0)+(j-local_y/2.0)*(j-local_y/2.0)+(k-local_z/2.0)*(k-local_z/2.0));
                if(ir<=xprd/20.0){
                  xi[0][na0]=1.0;
                  xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
                }
                //write input file
                //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
                data_fftw[index].re = xi[0][na0];
                data_fftw[index].im = xi[0][na1];
              }
              else if(((k-(local_z/4.0))<eps) &&                              // {100} plane
              ((k-(local_z/4.0))>-eps) &&
              (is==0)){
                ir = sqrt(((lxs+i)-N1/2.0)*((lxs+i)-N1/2.0)+(j-local_y/2.0)*(j-local_y/2.0)+(k-local_z/4.0)*(k-local_z/4.0));

                if(ir<=xprd/20.0){
                  xi[0][na0]=1.0;
                  xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
                }
                //write input file
                //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
                data_fftw[index].re = xi[0][na0];
                data_fftw[index].im = xi[0][na1];
              }
            }
          }
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }

        /* -----------------------------------------------------------------------
        Initial Crystal mode1
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::initial_sxtal_StraightOrtho()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          double xf=0, yf=0,zf=0;
          double eps=2e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
              for(j=0;j<N2;j++)
        	for(k=0;k<N3;k++){
        	  na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
        	  na = 2*(i*N2*N3+j*N3+k);
        	  index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
        	  na1 = na0+1;
        	  xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
        	  xi[0][na1] = 0.0; // just initializing them to 0.. not asking for an input
        	  xo[index] = 0.0;

            if(((k-local_z/2.0)<eps) &&     // {001} plane
              ((k-local_z/2.0)>-eps) &&
              ((lxs+i)>=(N1*0.25)) && ((lxs+i)<=(N1*3.0/4.0)) &&  //(j>(N2/4.0)) && (j<(3.0*N2/4.0)) && ////
              (is==0)){


              xi[0][na0]=1.0;
              xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
            }
              //write input file
              //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
              data_fftw[index].re = xi[0][na0];
              data_fftw[index].im = xi[0][na1];

        	}
          }
          if (app->nphase == 2) set_phase_2Px();
          if (app->nphase == 3) set_phase_3Px();
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }
        /* -----------------------------------------------------------------------
        set 2 phases
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::set_phase_2Px()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          double xf=0, yf=0,zf=0;
          double eps=2e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing phases ...\n");
            if (screen) fprintf(screen,"Initializing phases ...\n");
          }

          is=0;

            for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
            index = i*N2*N3+j*N3+k;
            phase[index] = 0;
           // sphere or cylin
            //ir = sqrt(((lxs+i)-N1/2.0)*((lxs+i)-N1/2.0)+(j-local_y/2.0)*(j-local_y/2.0)+(k-local_z/2.0)*(k-local_z/2.0));
            //ir = sqrt(((lxs+i)-N1/2.0)*((lxs+i)-N1/2.0)+(j-local_y/2.0)*(j-local_y/2.0));
	          /*ir = sqrt(((lxs+i)-N1/2.0)*((lxs+i)-N1/2.0)+(k-local_z/2.0)*(k-local_z/2.0));
             if(ir<=16.0){
               phase[index] = 1;
             }*/

            // slab
            //if(((lxs+i)>=28.0) && ((lxs+i)<=36.0) && (k>=8) && (k<=24)){
            if(((lxs+i)>=28.0) && ((lxs+i)<=36.0) && (j>=8) && (j<=23) && (k>=8) && (k<=23)){
            //if(((lxs+i)>=32.0) ){
              phase[index] = 1;
            } 

          }

          for(im=0;im<dimension;im++) {
            for(jm=0;jm<dimension;jm++){
              for(i=0;i<lN1;i++) {
                for(j=0;j<N2;j++) {
                  for(k=0;k<N3;k++) {
                    na = i*N2*N3+j*N3+k;
                    na0 = 2*(i*N2*N3+j*N3+k+im*lN1*N2*N3+jm*lN1*N2*N3*dimension);
                    na1 = na0+1;
                    index = i*N2*N3+j*N3+k+im*lN1*N2*N3+jm*lN1*N2*N3*dimension;
                    if((phase[na] == 1)){//&&((k-local_z/2.0)<eps)&&((k-local_z/2.0)>-eps)){
                      vstrain[na0] = 0.001;
                    }
                    vstrain_complex[index].re = vstrain[na0];
                    vstrain_complex[index].im = vstrain[na1];
                  }
                }
              }
            }
          }

          if(me==0){
            if (logfile) fprintf(logfile,"set_phase_2Px COMPLETE ...\n");
            if (screen) fprintf(screen,"set_phase_2Px COMPLETE  ...\n");
          }
          return;
        }
        /* -----------------------------------------------------------------------
        set 2 phases
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::set_phase_3Px()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          double xf=0, yf=0,zf=0;
          double eps=2e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing phases ...\n");
            if (screen) fprintf(screen,"Initializing phases ...\n");
          }

          is=0;

            for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
            index = i*N2*N3+j*N3+k;
            phase[index] = 0;
            if(((lxs+i)>=4.0) && ((lxs+i)<=27.0) && (j>=8) && (j<=23) && (k>=8) && (k<=23)){
              phase[index] = 1;
            }
            if(((lxs+i)>=36.0) && ((lxs+i)<=59.0) && (j>=8) && (j<=23) && (k>=8) && (k<=23)){
              phase[index] = 2;
            }
          }

          for(im=0;im<dimension;im++) {
            for(jm=0;jm<dimension;jm++){
              for(i=0;i<lN1;i++) {
                for(j=0;j<N2;j++) {
                  for(k=0;k<N3;k++) {
                    na = i*N2*N3+j*N3+k;
                    na0 = 2*(i*N2*N3+j*N3+k+im*lN1*N2*N3+jm*lN1*N2*N3*dimension);
                    na1 = na0+1;
                    index = i*N2*N3+j*N3+k+im*lN1*N2*N3+jm*lN1*N2*N3*dimension;
                    if((phase[na] != 0)){//&&((k-local_z/2.0)<eps)&&((k-local_z/2.0)>-eps)){
                      vstrain[na0] = 0.001;
                    }
                    vstrain_complex[index].re = vstrain[na0];
                    vstrain_complex[index].im = vstrain[na1];
                  }
                }
              }
            }
          }

          if(me==0){
            if (logfile) fprintf(logfile,"set_phase_3Px COMPLETE ...\n");
            if (screen) fprintf(screen,"set_phase_3Px COMPLETE  ...\n");
          }
          return;
        }

        /* -----------------------------------------------------------------------
        Initial Crystal mode1
        ---------------------------------------------------------------------*/
        void FFTW_Slab_Conduct::initial_sxtal_StraightNonOrtho()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-1;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
              for(j=0;j<N2;j++)
        	for(k=0;k<N3;k++){
        	  na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
        	  na = 2*(i*N2*N3+j*N3+k);
        	  index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
        	  na1 = na0+1;
        	  xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
        	  xi[0][na1] = 0.0;
        	  xo[index] = 0.0;


        	  if(((k-local_z/2.0)<eps) &&     // {111} plane
        	     ((k-local_z/2.0)>-eps) &&
        	     (((2*(lxs+i)-j)<=(-N1)) || (((2*(lxs+i)-j)>=0.0) && ((2*(lxs+i)-j)<=(N1)))) &&
        	     (is==0)){

        	    xi[0][na0]=1.0;
        	    xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];

        	    data_fftw[index].re = xi[0][na0];
        	    data_fftw[index].im = xi[0][na1];
        	  }
        	}
          }
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }
        /* -----------------------------------------------------------------------
        Initial Crystal HCP Basal
        ---------------------------------------------------------------------*/
        //CLaire added 07/27/18 to pair with the new app_3d_hcp_basal.cpp
        void FFTW_Slab_Conduct::initial_sxtal_3D3ShcpBasal()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-2;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
              na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
              na = 2*(i*N2*N3+j*N3+k);
              index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
              na1 = na0+1;
              xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
              xi[0][na1] = 0.0;
              xo[index] = 0.0;

              // Rotate to grid
              // xf = ((lxs+i)*four[0][0] + j*four[0][1] + k*four[0][2])/normfour[0];
              // yf = ((lxs+i)*four[1][0] + j*four[1][1] + k*four[1][2])/normfour[1];
              // zf = ((lxs+i)*four[2][0] + j*four[2][1] + k*four[2][2])/normfour[2];

              // xf = ((lxs+i)*prim[0][0] + j*prim[0][1] + k*prim[0][2]);
              // yf = ((lxs+i)*prim[1][0] + j*prim[1][1] + k*prim[1][2]);
              // zf = ((lxs+i)*prim[2][0] + j*prim[2][1] + k*prim[2][2]);

              if(((k-local_z/2.0)<eps) &&     // {001} plane
              ((k-local_z/2.0)>-eps) &&
              ((lxs+i)>(N1/4.0)) && ((lxs+i)<(3.0*N1/4.0)) && //(j>(N2/4.0)) && (j<(3.0*N2/4.0))
              (is==0)){
                // if(((zf-cntf[2])<eps) &&     // {111} plane
                //    ((zf-cntf[2])>-eps) &&
                //    (xf<(cntf[0]+xprdf/6.0+eps)) && (xf>(cntf[0]-xprdf/6.0)) &&
                //    (is==0)){
                // if(((xf-cntf[0])+(yf-cntf[1])+(zf-cntf[2])<eps) &&     // {111} plane
                //    ((xf-cntf[0])+(yf-cntf[1])+(zf-cntf[2])>-eps) &&
                //    (is==1)){
                xi[0][na0]=1.0;
                xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];

                //write input file
                //fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[0][na0], xi[0][na1], xi_sum[0][na]);
                data_fftw[index].re = xi[0][na0];
                data_fftw[index].im = xi[0][na1];
              }
              // else if(((zf-(cntf[2]/2.0))<eps) &&                              // {100} plane
              // 	  ((zf-(cntf[2]/2.0))>-eps) &&
              // 	  (is==0)){
              //   ir = sqrt((xf-cntf[0])*(xf-cntf[0])+(yf-cntf[1])*(yf-cntf[1])+(zf-(cntf[2]/2.0))*(zf-(cntf[2]/2.0)));
              //   if(ir<=xprdf/10.0){
              //     xi[0][na0]=1.0;
              //     xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
              //   }
            } /* ijk */
          } /* is */
          //fclose(of0);
          if(me==0){
            if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
            if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
          }
          return;
        }

        /* -----------------------------------------------------------------------
        Initial Crystal HCP Notch
        ---------------------------------------------------------------------*/
        //CLaire added 09/09/18 to pair with the new app_3d3s_hcp_notch.cpp
        void FFTW_Slab_Conduct::initial_sxtal_3D3ShcpNotch()
        {
          int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
          int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
          double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
          FILE *of0;
          char infile[100], input[100], c[10];
          double xf=0, yf=0,zf=0;
          double eps=1e-2;
          double ir=0;

          int NS = slip_systems;
          int N1 = nx;
          int lN1 = local_x;
          int lxs = local_x_start;
          int N2 = local_y;
          int N3 = local_z;
          obsden = 0.1;

          if(me==0){
            if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
            if (screen) fprintf(screen,"Initializing order parameters ...\n");
          }
          for(is=0;is<NS;is++){
            for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
              na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
              na = 2*(i*N2*N3+j*N3+k);
              index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
              na1 = na0+1;
              xi[0][na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
              xi[0][na1] = 0.0;
              xo[index] = 0.0;

              if((lxs+i) < 5 || j < 5 || k < 5){
                //    if(((lxs+i)-5) < eps || (j-5) < eps || (k-5) < eps){
                //      if(is == 0 /*&& k >= 5*/){
                //        xi[0][na0] = 0.0;
                //      }
                xo[index] = 1.0;
              }
              if(((k-local_z/2.0)<eps) && ((k-local_z/2.0)>-eps)){    // {0001} plane
              //    if(k == (N3/2.0)){
              //      ir = (((lxs+i)-(N1-1))*((lxs+i)-(N1-1)))+((j-N2/2)*(j-N2/2));
              //      //if(ir<=N1/2){
              //      if(ir<=25){
              //        if((lxs+i)== (N1-1) || (lxs+i)== (N1-2)){
              //          xo[index] = 1.0;
              //          if(is == 0){
              //      xi[0][na0] = 1.0;
              //          }
              //        }
              //      }
              //ledge or notch
              ir = (((lxs+i)-5.0)*((lxs+i)-5.0))+((j-N2/2)*(j-N2/2));
              if(ir<=25){
                //if((lxs+i)>=5){
                if((lxs+i) == 5 || (lxs+i) == 6){
                  xo[index] = 1.0;
                  if(is == 0 /*&& k >= 5*/){
                    //if(is == 1 && k == (N3-1)/2 && j >= 7){
                    xi[0][na0] = 1.0;
                  }
                }
              } /* ir */
              //      xi[0][na0] = 1.0;
              //      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];

              data_fftw[index].re = xi[0][na0];
              data_fftw[index].im = xi[0][na1];
            }
          } //ijk
        } //is
        //fclose(of0);
        if(me==0){
          if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
          if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
        }
        //    MPI_Barrier(MPI_COMM_WORLD);
        return;
      }

      void FFTW_Slab_Conduct::initial_sxtal_frank_read()
      {
        int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
        int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
        double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
        FILE *of0;
        char infile[100], input[100], c[10];
        double xf=0, yf=0,zf=0;
        double eps=1e-1;

        int NS = slip_systems;
        int N1 = nx;
        int lN1 = local_x;
        int lxs = local_x_start;
        int N2 = local_y;
        int N3 = local_z;
        obsden = 0.1;

        int fr_halflen = 10;

        if(me==0){
          if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
          if (screen) fprintf(screen,"Initializing order parameters ...\n");
        }
        for(is=0;is<NS;is++)
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
          na = 2*(i*N2*N3+j*N3+k);
          index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
          na1 = na0+1;
          xi[0][na0] = 0.0;
          xi[0][na1] = 0.0;
          xo[index] = 0.0;
          data_fftw[index].re = 0.0;
          data_fftw[index].im = 0.0;

          /*if(is == 1 && (lxs+i) == N1/2 && k <= N3/2 && j < N2/2 + fr_halflen && j >= N2/2 - fr_halflen){
          xi[0][na0]=1.0;
          xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
        }
        if(k != N3/2){ //confine slip to mid z plane
        xo[index] = 1.0;
      }*/
      if(is == 1 && (lxs+i) < N1/2 + fr_halflen && (lxs+i) >= N1/2 - fr_halflen && j <= N2/2 && k == N3/2){
        xi[0][na0] = 1.0;
        xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
      }
      if(j != 0 && j != N2/2){
        xo[index] = 1.0;
      }

      data_fftw[index].re = xi[0][na0];
      data_fftw[index].im = xi[0][na1];
    }
    if(me==0){
      if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
      if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
    }
    return;
  }

  void FFTW_Slab_Conduct::initial_sxtal_3SBCC()
  {
    int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
    int na0, na, na1, index, indexm, ptindex, *nodes, rtn_val0, rtn_val1, layer;
    double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
    FILE *of0;
    char infile[100], input[100], c[10];
    double xf=0, yf=0,zf=0;
    double eps=1e-1;

    int NS = slip_systems;
    int N1 = nx;
    int lN1 = local_x;
    int lxs = local_x_start;
    int N2 = local_y;
    int N3 = local_z;
    obsden = 0.1;

    if(me==0){
      if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
      if (screen) fprintf(screen,"Initializing order parameters ...\n");
    }
    for(is=0;is<NS;is++)
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
      na = 2*(i*N2*N3+j*N3+k);
      index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
      ptindex = i*N2*N3+j*N3+k;
      na1 = na0+1;
      xi[0][na0] = 0.0;
      xi[0][na1] = 0.0;
      xo[index] = 0.0;
      data_fftw[index].re = 0.0;
      data_fftw[index].im = 0.0;

      if(is == 0 && (lxs+i) == 0 && k >= int(N3/4.0) && k <= int(N3*3.0/4.0)){
        xi[0][na0] = 1.0;
        xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
      }

      data_fftw[index].re = xi[0][na0];
      data_fftw[index].im = xi[0][na1];
    }
    if(me==0){
      if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
      if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
    }
    return;
  }
void FFTW_Slab_Conduct::initial_sxtal_3SFCC()
{
  int is, ip, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
  int na0, na, na1, ptindex, index, index1, index2, index3, indexm, *nodes, rtn_val0, rtn_val1, layer;
  double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
  FILE *of0;
  char infile[100], input[100], c[10];
  double xf=0, yf=0,zf=0;
  double eps=1e-1;

  int NS = slip_systems;
  int NP = num_planes;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  obsden = 0.1;
  int eff_vol_mpi = 0;

  if(me==0){
    if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
    if (screen) fprintf(screen,"Initializing order parameters ...\n");
  }
  for(is=0;is<NS;is++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na = 2*(i*N2*N3+j*N3+k);
    index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
    ptindex = i*N2*N3+j*N3+k;
    na1 = na0+1;
    xi[0][na0] = 0.0;
    xi[0][na1] = 0.0;
    xo[index] = 0.0;
    xcount[ptindex] = 0.0;
    data_fftw[index].re = 0.0;
    data_fftw[index].im = 0.0;

    if(is == 0 && (lxs+i) == 0 && k >= int(N3/4.0) && k <= int(N3*3.0/4.0)){
      xi[0][na0] = 1.0;
      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
    }
    data_fftw[index].re = xi[0][na0];
    data_fftw[index].im = xi[0][na1];
  }
  if(app->dflag == 1){
  for(ip=0;ip<NP;ip++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
    na1=na0+1;
    dreg[0][na0] = 0.0;
    dreg[0][na1] = 0.0;
    dreg[0][na0] = 2.0*xi[0][index1] - xi[0][index2] - xi[0][index3];
    dreg[0][na1] = 3.0*xi[0][index2] - 3.0*xi[0][index3];
    dregp[0][na0] = 0.0;
    dregp[0][na1] = 0.0;
    dregp[0][na0] = dreg[0][na0]*cos(M_PI/6.0)+dreg[0][na1]*sin(M_PI/6.0); // Trailing partial
    dregp[0][na1] = dreg[0][na0]*cos(M_PI/6.0)-dreg[0][na1]*sin(M_PI/6.0); // Leading partial
  }
  }
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    ptindex = i*N2*N3+j*N3+k;
    if (xcount[ptindex] > 0.0){
      eff_vol_mpi += 1;
    }
  }
  MPI_Allreduce(&eff_vol_mpi, &eff_vol, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(me==0){
    if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ... effective volume %d\n",eff_vol);
    if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ... effective volume %d\n",eff_vol);
  }
  return;
}
void FFTW_Slab_Conduct::initial_sxtal_3SFCCPlanar()
{
  int is, ip, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
  int na0, na, na1, ptindex, index, index1, index2, index3, indexm, *nodes, rtn_val0, rtn_val1, layer;
  double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
  FILE *of0;
  char infile[100], input[100], c[10];
  double xf=0, yf=0,zf=0;
  double eps=1e-1;

  int NS = slip_systems;
  int NP = num_planes;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  obsden = 0.1;
  int eff_vol_mpi = 0;

  if(me==0){
    if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
    if (screen) fprintf(screen,"Initializing order parameters ...\n");
  }
  for(is=0;is<NS;is++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na = 2*(i*N2*N3+j*N3+k);
    index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
    ptindex = i*N2*N3+j*N3+k;
    na1 = na0+1;
    xi[0][na0] = 0.0;
    xi[0][na1] = 0.0;
    xo[index] = 0.0;
    xcount[ptindex] = 0.0;
    data_fftw[index].re = 0.0;
    data_fftw[index].im = 0.0;

    if(is == 0 && (lxs+i) == 0 && k >= int(N3/4.0) && k <= int(N3*3.0/4.0)){
      xi[0][na0] = 1.0;
      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
    }
    if(is < 3 && (lxs+i) == 0){
      xo[index] = 0.0;
      xcount[ptindex] = 1.0;
    }
    else{
      xo[index] = 1.0;
    }
    data_fftw[index].re = xi[0][na0];
    data_fftw[index].im = xi[0][na1];
  }
  if(app->dflag == 1){
  for(ip=0;ip<NP;ip++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
    na1=na0+1;
    dreg[0][na0] = 0.0;
    dreg[0][na1] = 0.0;
    dreg[0][na0] = 2.0*xi[0][index1] - xi[0][index2] - xi[0][index3];
    dreg[0][na1] = 3.0*xi[0][index2] - 3.0*xi[0][index3];
    dregp[0][na0] = 0.0;
    dregp[0][na1] = 0.0;
    dregp[0][na0] = dreg[0][na0]*cos(M_PI/6.0)+dreg[0][na1]*sin(M_PI/6.0); // Trailing partial
    dregp[0][na1] = dreg[0][na0]*cos(M_PI/6.0)-dreg[0][na1]*sin(M_PI/6.0); // Leading partial
  }
  }
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    ptindex = i*N2*N3+j*N3+k;
    if (xcount[ptindex] > 0.0){
      eff_vol_mpi += 1;
    }
  }
  MPI_Allreduce(&eff_vol_mpi, &eff_vol, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(me==0){
    if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ... effective volume %d\n",eff_vol);
    if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ... effective volume %d\n",eff_vol);
  }
  return;
}
void FFTW_Slab_Conduct::initial_sxtal_6SFCCPlanar()
{
  int is, ip, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
  int na0, na, na1, ptindex, index, index1, index2, index3, indexm, *nodes, rtn_val0, rtn_val1, layer;
  double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
  FILE *of0;
  char infile[100], input[100], c[10];
  double xf=0, yf=0,zf=0;
  double eps=1e-1;

  int NS = slip_systems;
  int NP = num_planes;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  obsden = 0.1;

  if(me==0){
    if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
    if (screen) fprintf(screen,"Initializing order parameters ...\n");
  }
  for(is=0;is<NS;is++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na = 2*(i*N2*N3+j*N3+k);
    index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
    ptindex = i*N2*N3+j*N3+k;
    na1 = na0+1;
    xi[0][na0] = 0.0;
    xi[0][na1] = 0.0;
    xo[index] = 0.0;
    data_fftw[index].re = 0.0;
    data_fftw[index].im = 0.0;

    if(is == 0 && (lxs+i) == 0 && k >= int(N3/4.0) && k <= int(N3*3.0/4.0)){
      xi[0][na0] = 1.0;
      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
    }
    if(is < 3 && (lxs+i) == 0 && ((k >= int(N3/8.0)) && (k <= int(N3*(7.0/8.0))))){
      xo[index] = 0.0;
    }
    else if (is > 2 && ((k == int(N3/8.0)) || (k == int(N3*(7.0/8.0))))){
      xo[index] = 0.0;
    }
    else{
      xo[index] = 1.0;
    }
    data_fftw[index].re = xi[0][na0];
    data_fftw[index].im = xi[0][na1];
  }
  if(app->dflag == 1){
  for(ip=0;ip<NP;ip++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
    na1=na0+1;
    dreg[0][na0] = 0.0;
    dreg[0][na1] = 0.0;
    dreg[0][na0] = 2.0*xi[0][index1] - xi[0][index2] - xi[0][index3];
    dreg[0][na1] = 3.0*xi[0][index2] - 3.0*xi[0][index3];
    dregp[0][na0] = 0.0;
    dregp[0][na1] = 0.0;
    dregp[0][na0] = dreg[0][na0]*cos(M_PI/6.0)+dreg[0][na1]*sin(M_PI/6.0); // Trailing partial
    dregp[0][na1] = dreg[0][na0]*cos(M_PI/6.0)-dreg[0][na1]*sin(M_PI/6.0); // Leading partial
  }
  }
  if(me==0){
    if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
    if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
  }
  return;
}
void FFTW_Slab_Conduct::initial_sxtal_6SFCCPlanar_Obs()
{
  int is, ip, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
  int na0, na, na1, ptindex, index, index1, index2, index3, indexm, *nodes, rtn_val0, rtn_val1, layer;
  double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
  FILE *of0;
  char infile[100], input[100], c[10];
  double xf=0, yf=0,zf=0;
  double eps=1e-1;

  int NS = slip_systems;
  int NP = num_planes;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  obsden = 0.1;
  int eff_vol_mpi = 0;

  int rank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Status status;

  if(me==0){
    if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
    if (screen) fprintf(screen,"Initializing order parameters ...\n");
  }
  for(is=0;is<NS;is++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na = 2*(i*N2*N3+j*N3+k);
    index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
    ptindex = i*N2*N3+j*N3+k;
    na1 = na0+1;
    xi[0][na0] = 0.0;
    xi[0][na1] = 0.0;
    xo[index] = 0.0;
    data_fftw[index].re = 0.0;
    data_fftw[index].im = 0.0;

    if(is == 0 && (lxs+i) == 0 && k >= int(N3/4.0) && k <= int(N3*3.0/4.0)){
      xi[0][na0] = 1.0;
      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
    }
    if(is < 3 && (lxs+i) == 0){
      xo[index] = 0.0;
      xcount[ptindex] = 1.0;
    }
    else if (is > 2 && ((k == int(N3/8.0)) || (k == int(N3*(7.0/8.0))))){
      xo[index] = 0.0;
      xcount[ptindex] = 1.0;
    }
    else{
      xo[index] = 1.0;
    }
    data_fftw[index].re = xi[0][na0];
    data_fftw[index].im = xi[0][na1];
  }
  if(app->dflag == 1){
  for(ip=0;ip<NP;ip++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
    na1=na0+1;
    dreg[0][na0] = 0.0;
    dreg[0][na1] = 0.0;
    dreg[0][na0] = 2.0*xi[0][index1] - xi[0][index2] - xi[0][index3];
    dreg[0][na1] = 3.0*xi[0][index2] - 3.0*xi[0][index3];
    dregp[0][na0] = 0.0;
    dregp[0][na1] = 0.0;
    dregp[0][na0] = dreg[0][na0]*cos(M_PI/6.0)+dreg[0][na1]*sin(M_PI/6.0); // Trailing partial
    dregp[0][na1] = dreg[0][na0]*cos(M_PI/6.0)-dreg[0][na1]*sin(M_PI/6.0); // Leading partial
  }
  }
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    ptindex = i*N2*N3+j*N3+k;
    if (xcount[ptindex] > 0.0){
      eff_vol_mpi += 1;
    }
  }
  MPI_Allreduce(&eff_vol_mpi, &eff_vol, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(me==0){
    if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ... effective volume %d\n",eff_vol);
    if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ... effective volume %d\n",eff_vol);
  }
  return;
}
void FFTW_Slab_Conduct::initial_sxtal_6SFCCPlanar_ObsInf()
{
  int is, ip, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
  int na0, na, na1, ptindex, index, index1, index2, index3, indexm, *nodes, rtn_val0, rtn_val1, layer;
  double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
  FILE *of0;
  char infile[100], input[100], c[10];
  double xf=0, yf=0,zf=0;
  double eps=1e-1;

  int NS = slip_systems;
  int NP = num_planes;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  obsden = 0.1;
  int eff_vol_mpi = 0;

  int rank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Status status;

  if(me==0){
    if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
    if (screen) fprintf(screen,"Initializing order parameters ...\n");
  }
  for(is=0;is<NS;is++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na = 2*(i*N2*N3+j*N3+k);
    index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
    ptindex = i*N2*N3+j*N3+k;
    na1 = na0+1;
    xi[0][na0] = 0.0;
    xi[0][na1] = 0.0;
    xo[index] = 0.0;
    data_fftw[index].re = 0.0;
    data_fftw[index].im = 0.0;

    if(is == 0 && (lxs+i) == 0 && k >= int(N3/4.0) && k <= int(N3*3.0/4.0)){
      xi[0][na0] = 1.0;
      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
    }
    if(is < 3 && (lxs+i) == 0 && ((k >= int(N3/8.0)) && (k <= int(N3*(7.0/8.0))))){
      xo[index] = 0.0;
      xcount[ptindex] = 1.0;
    }
    else if (is > 2 && ((k == int(N3/8.0)) || (k == int(N3*(7.0/8.0))))){
      xo[index] = 0.0;
      xcount[ptindex] = 1.0;
    }
    else{
      xo[index] = 1.0;
    }
    data_fftw[index].re = xi[0][na0];
    data_fftw[index].im = xi[0][na1];
  }
  if(app->dflag == 1){
  for(ip=0;ip<NP;ip++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
    na1=na0+1;
    dreg[0][na0] = 0.0;
    dreg[0][na1] = 0.0;
    dreg[0][na0] = 2.0*xi[0][index1] - xi[0][index2] - xi[0][index3];
    dreg[0][na1] = 3.0*xi[0][index2] - 3.0*xi[0][index3];
    dregp[0][na0] = 0.0;
    dregp[0][na1] = 0.0;
    dregp[0][na0] = dreg[0][na0]*cos(M_PI/6.0)+dreg[0][na1]*sin(M_PI/6.0); // Trailing partial
    dregp[0][na1] = dreg[0][na0]*cos(M_PI/6.0)-dreg[0][na1]*sin(M_PI/6.0); // Leading partial
  }
  }
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    ptindex = i*N2*N3+j*N3+k;
    if (xcount[ptindex] > 0.0){
      eff_vol_mpi += 1;
    }
  }
  MPI_Allreduce(&eff_vol_mpi, &eff_vol, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(me==0){
    if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ... effective volume %d\n",eff_vol);
    if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ... effective volume %d\n",eff_vol);
  }
  return;
}
void FFTW_Slab_Conduct::initial_sxtal_6SFCCPlanarKink()
{
  int is, ip, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
  int na0, na, na1, index, index1, index2, index3, indexm, *nodes, rtn_val0, rtn_val1, layer;
  double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
  FILE *of0;
  char infile[100], input[100], c[10];
  double xf=0, yf=0,zf=0;
  double eps=1e-1;

  int NS = slip_systems;
  int NP = num_planes;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  obsden = 0.1;

  if(me==0){
    if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
    if (screen) fprintf(screen,"Initializing order parameters ...\n");
  }
  for(is=0;is<NS;is++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na = 2*(i*N2*N3+j*N3+k);
    index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
    na1 = na0+1;
    xi[0][na0] = 0.0;
    xi[0][na1] = 0.0;
    xo[index] = 0.0;
    data_fftw[index].re = 0.0;
    data_fftw[index].im = 0.0;

    if(is == 0 && (lxs+i) == 0 && k >= int(N3/4.0) && k <= int(N3*3.0/4.0)){
      xi[0][na0] = 1.0;
      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
    }
    if(is == 0 && (lxs+i) == 0 && (j >= int(N2/2.0)-5 && j <= int(N2/2.0)+5) && (k == int(N3/4.0)-1 || k == int(N3*3.0/4.0)+1)){
      xi[0][na0] = 1.0;
      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
    }
    if(is < 3 && (lxs+i) == 0){
      xo[index] = 0.0;
    }
    else if (is > 2 && ((k == int(N3/4.0)) || (k == int(N3*3.0/4.0)))){
      xo[index] = 0.0;
    }
    else{
      xo[index] = 1.0;
    }

    data_fftw[index].re = xi[0][na0];
    data_fftw[index].im = xi[0][na1];
  }
  if(app->dflag == 1){
  for(ip=0;ip<NP;ip++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
    na1=na0+1;
    dreg[0][na0] = 0.0;
    dreg[0][na1] = 0.0;
    dreg[0][na0] = 2.0*xi[0][index1] - xi[0][index2] - xi[0][index3];
  }
  }
  if(me==0){
    if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
    if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
  }
  return;
}
void FFTW_Slab_Conduct::initial_sxtal_6SFCCGPZone()
{
  int is, ip, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
  int na0, na, na1, ptindex, index, index1, index2, index3, indexm, *nodes, rtn_val0, rtn_val1, layer;
  double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
  FILE *of0;
  char infile[100], input[100], c[10];
  double xf=0, yf=0,zf=0;
  double eps=1e-1;

  int NS = slip_systems;
  int NP = num_planes;
  int N1 = nx;
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  obsden = 0.1;

  if(me==0){
    if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
    if (screen) fprintf(screen,"Initializing order parameters ...\n");
  }
  for(is=0;is<NS;is++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na = 2*(i*N2*N3+j*N3+k);
    index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
    ptindex = i*N2*N3+j*N3+k;
    na1 = na0+1;
    xi[0][na0] = 0.0;
    xi[0][na1] = 0.0;
    xo[index] = 0.0;
    data_fftw[index].re = 0.0;
    data_fftw[index].im = 0.0;

    if(is == 0 && (lxs+i) == int(N1/2)-1 && k >= int(N3/4.0) && k <= int(N3*3.0/4.0)){
      xi[0][na0] = 1.0;
      xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
    }
    // if(is < 3 && (lxs+i) == int(N1/2)-1/* && ((k >= int(N3/8.0)) && (k <= int(N3*7.0/8.0)))*/){
    //   xo[index] = 0.0;
    // }
    // else if (is > 2 && ((k == int(N3/8.0)) || (k == int(N3*7.0/8.0)))){
    //   xo[index] = 0.0;
    // }
    // else{
    //   xo[index] = 1.0;
    // }
    int gprd = 4;
    int gpct = 8;

    // for(int gp=0;gp<gpct*2;gp++){
    //   if(gp % 2 != 0){
    //     if(((lxs+i+k) == int(N1/2)-1+int(N3/8) || (lxs+i+k) == int(N1/2)-1+int(N3*7/8)) && (lxs+i-int(N1/2))*(lxs+i-int(N1/2)) + (j-int(N2*gp/(gpct*2)))*(j-int(N2*gp/(gpct*2))) <= gprd*gprd) xo[index] = 2.0; // (001)
    //   }
    // }

    for(int gp=0;gp<gpct*2;gp++){
      if(gp % 2 != 0){
        if(((j+k) == int(N2*gp/(gpct*2))+int(N3/8) || (j+k) == int(N2*gp/(gpct*2))+int(N3*7/8)) && (lxs+i-int(N1/2))*(lxs+i-int(N1/2)) + (j-int(N2*gp/(gpct*2)))*(j-int(N2*gp/(gpct*2))) <= gprd*gprd) xo[index] = 2.0; // (010)
      }
    }

    // for(int gp=0;gp<gpct*2;gp++){
    //   if(gp % 2 != 0){
    //     if(((lxs+i+j) == int(N1/2)-1+int(N2*gp/(gpct*2))) && (lxs+i-int(N1/2)+1)*(lxs+i-int(N1/2)+1) + (j-int(N2*gp/(gpct*2)))*(j-int(N2*gp/(gpct*2))) + (k-int(N3/8.0))*(k-int(N3/8.0)) <= gprd*gprd) xo[index] = 2.0; // (100)
    //     else if(((lxs+i+j) == int(N1/2)-1+int(N2*gp/(gpct*2))) && (lxs+i-int(N1/2)+1)*(lxs+i-int(N1/2)+1) + (j-int(N2*gp/(gpct*2)))*(j-int(N2*gp/(gpct*2))) + (k-int(N3*7.0/8.0))*(k-int(N3*7.0/8.0)) <= gprd*gprd) xo[index] = 2.0; // (100)
    //   }
    // }

    // opc[is][ptindex] = std::complex<double> (xi[0][na0],xi[0][na1]);
    data_fftw[index].re = xi[0][na0];
    data_fftw[index].im = xi[0][na1];
  }
  if(app->dflag == 1){
  for(ip=0;ip<NP;ip++)
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1=2*(i*N2*N3+j*N3+k+(3*ip+0)*lN1*N2*N3);
    index2=2*(i*N2*N3+j*N3+k+(3*ip+1)*lN1*N2*N3);
    index3=2*(i*N2*N3+j*N3+k+(3*ip+2)*lN1*N2*N3);
    na0 = 2*(i*N2*N3+j*N3+k+ip*lN1*N2*N3);
    na1=na0+1;
    dreg[0][na0] = 0.0;
    dreg[0][na1] = 0.0;
    dreg[0][na0] = 2.0*xi[0][index1] - xi[0][index2] - xi[0][index3];
    dreg[0][na1] = 3.0*xi[0][index2] - 3.0*xi[0][index3];
    dregp[0][na0] = 0.0;
    dregp[0][na1] = 0.0;
    dregp[0][na0] = dreg[0][na0]*cos(M_PI/6.0)+dreg[0][na1]*sin(M_PI/6.0); // Trailing partial
    dregp[0][na1] = dreg[0][na0]*cos(M_PI/6.0)-dreg[0][na1]*sin(M_PI/6.0); // Leading partial
  }
  }
  if(me==0){
    if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
    if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
  }
  return;
}
  void FFTW_Slab_Conduct::initial_sxtal_loop()
  {
    int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
    int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
    double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
    FILE *of0;
    char infile[100], input[100], c[10];
    double xf=0, yf=0,zf=0;
    double eps=1e-1;

    int NS = slip_systems;
    int N1 = nx;
    int lN1 = local_x;
    int lxs = local_x_start;
    int N2 = local_y;
    int N3 = local_z;
    int dist[3];
    obsden = 0.1;

    if(me==0){
      if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
      if (screen) fprintf(screen,"Initializing order parameters ...\n");
    }
    for(is=0;is<NS;is++)
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
      na = 2*(i*N2*N3+j*N3+k);
      index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
      na1 = na0+1;
      xi[0][na0] = 0.0;
      xi[0][na1] = 0.0;
      xo[index] = 0.0;
      data_fftw[index].re = 0.0;
      data_fftw[index].im = 0.0;
      dist[0] = -1*(j-N2/2)+(k-N3/2);
      dist[1] = (j-N2/2)-(k-N3/2);
      dist[2] = (j-N2/2)+(k-N3/2);

      if(is == 0 && (lxs+i) == 0 && dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2] < 3*16*16){
        xi[0][na0] = 1.0;
        xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
      }

      data_fftw[index].re = xi[0][na0];
      data_fftw[index].im = xi[0][na1];
    }
    if(me==0){
      if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
      if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
    }
    return;
  }

  void FFTW_Slab_Conduct::initial_sxtal_cubic_fcc_edge()
  {
    int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
    int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
    double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
    FILE *of0;
    char infile[100], input[100], c[10];
    double xf=0, yf=0,zf=0;
    double eps=1e-1;

    int NS = slip_systems;
    int N1 = nx;
    int lN1 = local_x;
    int lxs = local_x_start;
    int N2 = local_y;
    int N3 = local_z;
    obsden = 0.1;

    if(me==0){
      if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
      if (screen) fprintf(screen,"Initializing order parameters ...\n");
    }
    for(is=0;is<NS;is++)
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
      na = 2*(i*N2*N3+j*N3+k);
      index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
      na1 = na0+1;
      xi[0][na0] = 0.0;
      xi[0][na1] = 0.0;
      data_fftw[index].re = 0.0;
      data_fftw[index].im = 0.0;
      xo[index] = 0.0;
      // In order to make the slip plane periodic we need two slip planes in the same domain that must be cubic
      if( is == 0 && ( ( ((lxs+i + j + k) >= ((N1-1)*0.5 + (N2-1)*0.5 + N3 - eps)) && ((lxs+i + j + k) <= ((N1-1)*0.5 + (N2-1)*0.5 + N3 + eps)) ) ||
          ( ((lxs+i + j + k) >= ((N1-1)*0.5 + (N2-1)*0.5 - eps)) && ((lxs+i + j + k) <= ((N1-1)*0.5 + (N2-1)*0.5 + eps)) ) )
          && ( ((lxs+i - j) >= ((N2+2)*3/4 - eps)) || ((lxs+i - j) <= (-(N2+2)*3/4 + eps))
              || (((lxs+i - j) > (-(N2)*0.25 + eps )) && ((lxs+i - j) < ((N2)*0.25 - eps ))) ) ){
                xi[0][na0] = 1.0;
                xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
      }
      data_fftw[index].re = xi[0][na0];
      data_fftw[index].im = xi[0][na1];
    }
    if(me==0){
      if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
      if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
    }
    return;
  }

  void FFTW_Slab_Conduct::initial_sxtal_cubic_fcc_screw()
  {
    int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
    int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
    double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
    FILE *of0;
    char infile[100], input[100], c[10];
    double xf=0, yf=0,zf=0;
    double eps=1e-1;

    int NS = slip_systems;
    int N1 = nx;
    int lN1 = local_x;
    int lxs = local_x_start;
    int N2 = local_y;
    int N3 = local_z;
    obsden = 0.1;

    if(me==0){
      if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
      if (screen) fprintf(screen,"Initializing order parameters ...\n");
    }
    for(is=0;is<NS;is++)
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
      na = 2*(i*N2*N3+j*N3+k);
      index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
      na1 = na0+1;
      xi[0][na0] = 0.0;
      xi[0][na1] = 0.0;
      data_fftw[index].re = 0.0;
      data_fftw[index].im = 0.0;
      xo[index] = 0.0;
      // In order to make the slip plane periodic we need two slip planes in the same domain that must be cubic
      if( is == 0 && ( ( ((lxs+i + j + k) >= ((N1-1)*0.5 + (N2-1)*0.5 + N3 - eps)) && ((lxs+i + j + k) <= ((N1-1)*0.5 + (N2-1)*0.5 + N3 + eps)) ) ||
          ( ((lxs+i + j + k) >= ((N1-1)*0.5 + (N2-1)*0.5 - eps)) && ((lxs+i + j + k) <= ((N1-1)*0.5 + (N2-1)*0.5 + eps)) ) )
          && ( ( k >= ( N3*0.25 - eps ) ) && ( k <= ( N2*3.0/4.0 + eps ) ) ) ){
              xi[0][na0] = 1.0;
              xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
      }
      data_fftw[index].re = xi[0][na0];
      data_fftw[index].im = xi[0][na1];
    }
    if(me==0){
      if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
      if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
    }
    return;
  }

  void FFTW_Slab_Conduct::initial_sxtal_edge()
  {
    int is, i, j, k, im, jm, km, ism, ia, ib, indexgb, grainb, num, count;
    int na0, na, na1, index, indexm, *nodes, rtn_val0, rtn_val1, layer;
    double nlx, nly, c0, c1, a, alpha, zeta, eta, d;
    FILE *of0;
    char infile[100], input[100], c[10];
    double xf=0, yf=0,zf=0;
    double eps=1e-1;

    int NS = slip_systems;
    int N1 = nx;
    int lN1 = local_x;
    int lxs = local_x_start;
    int N2 = local_y;
    int N3 = local_z;
    obsden = 0.1;

    if(me==0){
      if (logfile) fprintf(logfile,"Initializing order parameters ...\n");
      if (screen) fprintf(screen,"Initializing order parameters ...\n");
    }
    for(is=0;is<NS;is++)
    for(i=0;i<lN1;i++)
    for(j=0;j<N2;j++)
    for(k=0;k<N3;k++){
      na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
      na = 2*(i*N2*N3+j*N3+k);
      index = i*N2*N3+j*N3+k+is*lN1*N2*N3;
      na1 = na0+1;
      xi[0][na0] = 0.0;
      xi[0][na1] = 0.0;
      xo[index] = 0.0;
      //N3 must be 3*N2 for this to make an edge dislocation !!!
      if(is == 0 && (lxs+i) == 0 && ((k <= 3*j && (k+N3/2 > 3*j)) || (k - N3/2) > 3*j)){
        xi[0][na0] = 1.0;
        xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];
      }
      data_fftw[index].re = xi[0][na0];
      data_fftw[index].im = xi[0][na1];
    }
    if(me==0){
      if (logfile) fprintf(logfile,"Initial configuration of order parameters COMPLETE ...\n");
      if (screen) fprintf(screen,"Initial configuration of order parameters COMPLETE  ...\n");
    }
    return;
  }

  /* -----------------------------------------------------------------------
  average strain
  ---------------------------------------------------------------------*/
  void FFTW_Slab_Conduct::stressfree_strain()
  {
    int ND = dimension;
    int N1 = nx;
    double size = static_cast<double>(N1);
    int lN1 = local_x;
    int lxs = local_x_start;
    int N2 = local_y;
    int N3 = local_z;
    int NP = num_planes;
    int NS = slip_systems;
    int nsize = N1*N2*N3;
    int i, j, k, l, is, nb, na0, ida, idb;
    // double S[ND][ND][ND][ND];
    double S44 = material->S44;
    double S12 = material->S12;
    double mu = material->mu;

    /*set S matrix*/
    // for (i=0; i<ND; i++)
    // for (j=0; j<ND; j++)
    // for (k=0; k<ND; k++)
    // for (l=0; l<ND; l++)
    // {
    //   S[i][j][k][l] = S44/4 * (DELTA(i,k)*DELTA(j,l)+DELTA(i,l)*DELTA(j,k))+S12*DELTA(i,j)*DELTA(k,l);
    // }


    /*calculating average stress free strain, avepsd*/
    for(is=0;is<NS;is++){
      for (ida=0; ida<ND; ida++){
        for (idb=0; idb<ND; idb++){
          for(i=0;i<lN1;i++){
            for(j=0;j<N2;j++){
              for(k=0;k<N3;k++){
                nb = 2*(k + j*N3 + i*N3*N2+ is*lN1*N2*N3);
                avepsd[ida][idb]  += eps[is][ida][idb] * xi[0][nb];  // integral 1/V Int_ \epsilon^0
              }
              //printf("avepsd[ida][idb] %lf, me %d, ida %d, idb %d\n", avepsd[ida][idb], me, ida, idb);
            }
          }

          MPI_Allreduce(&avepsd[ida][idb], &ave_epsd[ida][idb], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          //printf("ave_epsd[ida][idb] %lf, me %d, ida %d, idb %d, ave_eps/nisze %lf\n", ave_epsd[ida][idb], me, ida, idb, ave_epsd[ida][idb]/nsize);
          ave_epsd[ida][idb] = ave_epsd[ida][idb]/nsize;
        }
      }
    }
    if (app->nphase>1) {
      for(ida=0;ida<ND;ida++){
        for (idb=0;idb<ND;idb++){
          for(i=0;i<lN1;i++){
            for(j=0;j<N2;j++){
              for (k=0;k<N3;k++){
                na0 = 2*(k + j*N3 + i*N2*N3 + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
                aveVs[ida][idb] += vstrain[na0]; //aveps only appears here to get total average strain in the box
              }
            }
          }
          MPI_Allreduce(&aveVs[ida][idb], &ave_Vs[ida][idb], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          ave_Vs[ida][idb] /= ((double) nsize);
        } //idb
      } //ida
    }
    /*calculating microscopic strain, avepst*/

    for(i=0;i<ND;i++){
      for(j=0;j<ND;j++){
        for(k=0;k<ND;k++){
          for(l=0;l<ND;l++){
            avepst[i][j] += S[i][j][k][l]*sigma_rot[k][l]*mu;    // Compliance & stress
            avepsts[i][j] += S[i][j][k][l]*sigma_rot[k][l]*mu;
            //avepst[i][j] += S[i][j][k][l]*sigma[k][l]*mu;    // Compliance & stress
            //avepsts[i][j] += S[i][j][k][l]*sigma[k][l]*mu;
            //printf("avepst[i][j] %10.0lf, me %d, i %d, i %d, S[i][j][k][l] %10.0lf, k %d, l %d, mu %lf, sigma[k][l] %lf\n", avepst[i][j], me, i, j, S[i][j][k][l], k, l, mu, sigma[k][l]);
            //avepst[i][j] += S[i][j][k][l]*sigma[k][l];
          }
        }
        avepst[i][j] += ave_epsd[i][j];   //epsilon bar Wang 2002
        if (app->nphase>1) {avepst[i][j] += ave_Vs[i][j];}
      }
      /*
      if(me==0){
      if (logfile) fprintf(logfile, "avepst[%d][j] = %lf, %lf, %lf\n",i, avepst[i][0],avepst[i][1],avepst[i][2]);
      if (screen) fprintf(screen, "avepst[%d][j] = %lf, %lf, %lf\n",i, avepst[i][0],avepst[i][1],avepst[i][2]);
    }
    */
  }

  return;
}

/* -----------------------------------------------------------------------
average strain
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::total_average_strain()
{
  int ND = dimension;
  int N1 = nx;
  double size = static_cast<double>(N1);
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  int nsize = N1*N2*N3;
  int i, j, k, l, is, nb, ida, idb, na0;

  for(ida=0;ida<ND;ida++){
    for (idb=0;idb<ND;idb++){
      aveps[ida][idb] = 0.0;
      aveps_NoExt[ida][idb] = 0.0;
    }
  }
  for(ida=0;ida<ND;ida++){
    for (idb=0;idb<ND;idb++){
      for(i=0;i<lN1;i++){
        for(j=0;j<N2;j++){
          for (k=0;k<N3;k++){
            na0 = 2*(k + j*N3 + i*N2*N3 + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
            aveps[ida][idb] += data_eps[na0]; //aveps only appears here to get total average strain in the box
            aveps_NoExt[ida][idb] += data_eps_NoExt[na0];
          }
        }
      }
      MPI_Reduce(&aveps[ida][idb], &ave_eps[ida][idb], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      ave_eps[ida][idb] /= nsize;
      MPI_Reduce(&aveps_NoExt[ida][idb], &ave_eps_NoExt[ida][idb], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      ave_eps_NoExt[ida][idb] /= nsize;
    } //idb
  } //ida

}

/* -----------------------------------------------------------------------
strain
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::strain()
{
  int ND = dimension;
  int N1 = nx;
  double size = static_cast<double>(N1);
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  int nsize = N1*N2*N3;
  int i, j, k, l, is, na0, na1, nb, psys, ida, idb, idc, idd, index, index2;
  int na11, na12, na13, na21, na22, na23, na31, na32, na33, ia, ib;
  fftwnd_mpi_plan iiplan;

  iiplan = fftw3d_mpi_create_plan(world, N1, N2, N3, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (i=0; i<lN1*N2*N3*ND*ND; i++){
    data_strain[i].re = 0.0;
    data_strain[i].im = 0.0;
    data_strain_NoExt[i].re = 0.0;
    data_strain_NoExt[i].im = 0.0;
  }
  for (i=0; i<lN1*N2*N3; i++){
    work_strain[i].re = 0.0;
    work_strain[i].im = 0.0;
  }
  for (i=0; i<2*lN1*N2*N3*ND*ND; i++){
    data_eps[i] = 0.0;
    data_eps_NoExt[i] = 0.0;
  }

  /*calculate the total strain */

  for(is=0;is<NS;is++){
    for (ida=0; ida<ND; ida++){
      for (idb=0; idb<ND; idb++){
        for(i=0;i<lN1;i++){
          for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
              index = i*N2*N3 + j*N3 + k + is*lN1*N2*N3;
              index2 = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
              nb = k + j*N3 + i*N2*N3 + is*lN1*N2*N3 + ida*lN1*N2*N3*NS + idb*lN1*N2*N3*NS*ND;
              data_strain[index2].re += data_fftw[index].re * FF[nb];
              data_strain[index2].im += data_fftw[index].im * FF[nb];

            }
          }
        }
      }
    }
  }
  if (app->nphase>1) {
    for(i=0;i<lN1;i++){ // Multiply virtual strain by DD
      for(j=0;j<N2;j++){
        for(k=0;k<N3;k++) {
            for(ida=0;ida<ND;ida++){
              for(idb=0;idb<ND;idb++){
                for(idc=0;idc<ND;idc++){
                  for(idd=0;idd<ND;idd++){
                index2 = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
                nb  = i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND + idc*lN1*N2*N3*ND*ND + idd*lN1*N2*N3*ND*ND*ND;
                index  = i*N2*N3 + j*N3 + k + idc*lN1*N2*N3 + idd*lN1*N2*N3*ND;
                data_strain[index2].re += vstrain_complex[index].re * EE[nb];
                data_strain[index2].im += vstrain_complex[index].im * EE[nb];
              }
            }
              }
            }
        }
      }
    }
  }

  for(i=0;i<ND;i++){
    for (j=0;j<ND;j++){
      psys = i*lN1*N2*N3 + j*lN1*N2*N3*ND;
      fftwnd_mpi(iiplan, 1, data_strain+psys, work_strain, FFTW_NORMAL_ORDER); /* Inverse FFT (multiple)*/
    }
  }
  for (i=0; i<lN1*N2*N3*ND*ND; i++){
    data_strain[i].re = data_strain[i].re/(nsize);
    data_strain[i].im = data_strain[i].im/(nsize);
    data_strain_NoExt[i].re = data_strain[i].re;
    data_strain_NoExt[i].im = data_strain[i].im;
  }

  // calculate stress free strain
  stressfree_strain();
  //add in other two terms in strain (added already into avepst: epsilon bar 0 and S x sigma)
  // for(ida=0;ida<ND;ida++){
  //   for (idb=0;idb<ND;idb++){
  //     for(i=0;i<lN1;i++){
  //       for(j=0;j<N2;j++){
  //         for (k=0;k<N3;k++){
  //           index = k + j*N3 + i*N2*N3 + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
  //           data_strain[index].re += avepst[ida][idb];     // Adding epsilon bar Wang 2002
  //         }
  //       }
  //     }
  //   }
  // }

  for(ida=0;ida<ND;ida++){
    for (idb=0;idb<ND;idb++){
      for(i=0;i<lN1;i++){
        for(j=0;j<N2;j++){
          for (k=0;k<N3;k++){
            na0 = 2*(k + j*N3 + i*N2*N3 + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
            na1=na0+1;
            index = k + j*N3 + i*N2*N3 + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
            index2 = k + j*N3 + i*N2*N3 + idb*lN1*N2*N3 + ida*lN1*N2*N3*ND;
            //eps_ij = 1/2(dui/dxj + duj/dxi) --> ida=i, idb=j
            data_eps[na0] = (data_strain[index].re + data_strain[index2].re)/2.0;   // symmetric part of the strain
            data_eps[na1] = (data_strain[index].im + data_strain[index2].im)/2.0;

	          data_eps_NoExt[na0] = (data_strain_NoExt[index].re + data_strain_NoExt[index2].re)/2.0;   // symmetric part of the strain
	          data_eps_NoExt[na1] = (data_strain_NoExt[index].im + data_strain_NoExt[index2].im)/2.0;
          }
        }
      }
    }
  }

  fftwnd_mpi_destroy_plan(iiplan);

  return;
}

/* -----------------------------------------------------------------------
stress
---------------------------------------------------------------------*/
void FFTW_Slab_Conduct::stress()
{
  int ND = dimension;
  int N1 = nx;
  double size = static_cast<double>(N1);
  int lN1 = local_x;
  int lxs = local_x_start;
  int N2 = local_y;
  int N3 = local_z;
  int NP = num_planes;
  int NS = slip_systems;
  int nsize = N1*N2*N3;
  int i, j, k, l, m, ida, idb, na, nb, na0, is, ia, ib, index;
  int na11, na12, na13, na21, na22, na23, na31, na32, na33, layer;
  int tcount[ND][ND], ccount[ND][ND], t_count[ND][ND], c_count[ND][ND];
  //double C[ND][ND][ND][ND];
  double mu= material->mu;
  int ortho = 0;
  double dcs[ND][ND]; // Dislocation coordinate system

  // double C44 = material->C44;
  // double C12 = material->C12;
  // double C11 = material->C11;
  //
  // // Why do we have to redefine this?????
  //
  // mu = C44-(2.0*C44+C12-C11)/5.0;
  // ll = C12-(2.0*C44+C12-C11)/5.0;
  // young = mu*(3*ll+2*mu)/(ll+mu);
  // xnu = young/2.0/mu-1.0;

  for (i=0; i<2*lN1*N2*N3*ND*ND; i++){
    data_sigma[i] =0;
    data_sigma_NoExt[i] =0;
    data_asigma[i] =0;
    data_tsigma[i] =0;
    data_psigma[i] =0;
    data_usigma[i] =0;
    data_rsigma[i] = 0;
    data_epsd[i]=0;
  }

  // This is the same calculation as the stress free strain
  // We could do both calculations above

  //for(is=0;is<NSV;is++)
  for(is=0;is<NS;is++){
    for (ida=0; ida<ND; ida++){
      for (idb=0; idb<ND; idb++){
        for(i=0;i<lN1;i++){
          for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
              na = 2*(i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
              nb = 2*(k + j*N3 + i*N3*N2 + is*lN1*N2*N3);
              data_epsd[na] += eps[is][ida][idb] * xi[0][nb];
            }
          }
        }
      }
    }
  }
  for (ida=0;ida<ND;ida++){
    for (idb=0;idb<ND;idb++){
      tcount[ida][idb] = 0;
      ccount[ida][idb] = 0;
      for(i=0;i<lN1;i++){
        for(j=0;j<N2;j++){
          for(k=0;k<N3;k++){
            na = 2*(i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
            for (m=0;m<ND;m++){
              for (l=0;l<ND;l++){
                na0 = 2*(i*N2*N3 + j*N3 + k + m*lN1*N2*N3 + l*lN1*N2*N3*ND);
                index = i*N2*N3 + j*N3 + k + m*lN1*N2*N3 + l*lN1*N2*N3*ND;

                if (app->nphase>1) {
                  data_sigma[na] += C[ida][idb][m][l]*(data_strain[index].re - data_epsd[na0] - avepsts[m][l] - vstrain[na0]);
                } else {
                  data_sigma[na] += C[ida][idb][m][l]*(data_strain[index].re + avepst[m][l] - data_epsd[na0]);
                  data_sigma_NoExt[na] += C[ida][idb][m][l]*(data_strain[index].re - data_epsd[na0] - avepsts[m][l]);
                  data_asigma[na] += C[ida][idb][m][l]*(avepst[m][l]);
                  data_psigma[na] += -C[ida][idb][m][l]*data_epsd[na0];
                  data_tsigma[na] += C[ida][idb][m][l]*(data_strain[index].re - avepst[m][l]);
                }
                data_sigma[na+1] = 0.0;
                data_sigma_NoExt[na+1] = 0.0;
              }
            }
            //data_sigma[na]+=sigma[ida][idb]*mu;
            data_sigma[na]+=sigma_rot[ida][idb]*mu;
            local_sigma[ida][idb] += data_sigma[na];
            local_sigma_NoExt[ida][idb] += data_sigma_NoExt[na];
          }
        }
      }
      MPI_Allreduce(&local_sigma[ida][idb], &ave_sigma[ida][idb], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      ave_sigma[ida][idb] /= nsize;
      MPI_Allreduce(&local_sigma_NoExt[ida][idb], &ave_sigma_NoExt[ida][idb], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      ave_sigma_NoExt[ida][idb] /= nsize;
    }
  }

  if((fabs(fft->vol)>0.99) && (fabs(fft->vol)<1.01)) ortho =1;
  if (ortho){
    if(me==0){
      if (logfile) fprintf(logfile,"Orthogonal stress rotation ...\n");
      if (screen) fprintf(screen,"Orthogonal stress rotation  ...\n");
    }
    for(ida=0; ida<ND; ida++){
      for(idb=0; idb<ND; idb++){
        for(i=0;i<lN1;i++){
          for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
              na = 2*(i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
              for(m=0; m<ND; m++)
              for(l=0; l<ND; l++){
                na0 = 2*(i*N2*N3 + j*N3 + k + m*lN1*N2*N3 + l*lN1*N2*N3*ND);
                data_rsigma[na] += invprim[ida][m]*data_sigma[na0]*prim[l][idb];
              }
            }
          }
        }
      }
    }
  } else{
    if(me==0){
      if (logfile) fprintf(logfile,"Nonorthogonal stress rotation ...\n");
      if (screen) fprintf(screen,"Nonorthogonal stress rotation  ...\n");
    }
    // Here the rotation is defined by taking only one slip system
    dcs[0][0] = xb[0][1]*xn[0][2] - xb[0][2]*xn[0][1];
    dcs[0][1] = -(xb[0][0]*xn[0][2] - xb[0][2]*xn[0][0]);
    dcs[0][2] = xb[0][0]*xn[0][1] - xb[0][1]*xn[0][0];
    for(idb=0; idb<ND; idb++){
      dcs[1][idb] = xb[0][idb];
      dcs[2][idb] = xn[0][idb];
    }
    for(ida=0; ida<ND; ida++){
      for(idb=0; idb<ND; idb++){
        for(i=0;i<lN1;i++){
          for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
              na = 2*(i*N2*N3 + j*N3 + k + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
              for(m=0; m<ND; m++)
              for(l=0; l<ND; l++){
                na0 = 2*(i*N2*N3 + j*N3 + k + m*lN1*N2*N3 + l*lN1*N2*N3*ND);
                data_rsigma[na] += dcs[ida][m]*data_sigma[na0]*dcs[idb][l];
              }
            }
          }
        }
      }
    }
  }

}
