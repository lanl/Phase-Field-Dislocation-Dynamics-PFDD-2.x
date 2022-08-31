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

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "fft_fftw_slab.h"
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

FFTW_Slab::FFTW_Slab(PFDD_C *pfdd_p, int narg, char **arg) : FFT(pfdd_p,narg,arg)
{

  user_procgrid[0] = nprocs;
  user_procgrid[1] = user_procgrid[2] = 1;
}

/* ---------------------------------------------------------------------- */

FFTW_Slab::~FFTW_Slab()
{

}

/* ----------------------------------------------------------------------
initialize fftw
------------------------------------------------------------------------- */

void FFTW_Slab::init()
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
  }

  /* ----------------------------------------------------------------------
  setup fftw
  ------------------------------------------------------------------------- */

  void FFTW_Slab::setup()
  {
    rotate_stress();
    frec();
    stiffness();
    greens_function();
    Bmatrix();
    Fmatrix();
  }
  /* ----------------------------------------------------------------------
  setup global box
  assumes boxlo/hi are already set
  ------------------------------------------------------------------------- */

  void FFTW_Slab::create_plan()
  {

    /*create plan and iplan for fftw*/

    plan = fftw3d_mpi_create_plan(world, nx, ny, nz, FFTW_FORWARD, FFTW_ESTIMATE);

    iplan = fftw3d_mpi_create_plan(world, nx, ny, nz, FFTW_BACKWARD, FFTW_ESTIMATE);

  }

  /* ----------------------------------------------------------------------
  allocate
  ------------------------------------------------------------------------- */

  void FFTW_Slab::allocate()
  {

    slip_systems = app->slip_systems;
    num_planes = app->num_planes;
    dimension = app->dimension;
    memory->create(data_fftw,total_local_size*slip_systems,"data_fftw");
    memory->create(data_real,total_local_size*slip_systems,"data_real");
    memory->create(temp_data,total_local_size*slip_systems,"temp_data");
    memory->create(work,total_local_size,"work");
    memory->create(work_strain,total_local_size,"work_strain");
    memory->create(data_core,total_local_size*slip_systems,"data_core");
    memory->create(data_strain,total_local_size*dimension*dimension,"data_strain");

    memory->create(xi,norder,2.0*slip_systems*local_x*local_y*local_z,"xi");
    memory->create(xi_sum,norder,2.0*local_x*local_y*local_z,"xi_sum");
    for(int i=0; i<norder; i++){
      for(int j=0; j<2.0*slip_systems*local_x*local_y*local_z; j++)
      xi[i][j] = 0;
      for(int j=0; j<2.0*local_x*local_y*local_z; j++)
      xi_sum[i][j] = 0.0;
    }
    memory->create(xo,slip_systems*local_x*local_y*local_z,"xo");
    memory->create(fx,local_x*local_y*local_z,"fx");
    memory->create(fy,local_x*local_y*local_z,"fy");
    memory->create(fz,local_x*local_y*local_z,"fz");
    memory->create(f,slip_systems*local_x*local_y*local_z,"f");
    memory->create(r,dimension,"r");
    memory->create(C,dimension,dimension,dimension,dimension,"C");
    //memory->create(G,local_x,local_y,local_z,dimension,dimension,"G");

    memory->create(BB,slip_systems*slip_systems*local_x*local_y*local_z,"BB");
    memory->create(FF,slip_systems*total_local_size*dimension*dimension,"FF");
    memory->create(DD,slip_systems*total_local_size*dimension*dimension,"DD");

    //Grad allocation
    memory->create(gradx,num_planes*local_x*local_y*local_z,"gradx");
    memory->create(grady,num_planes*local_x*local_y*local_z,"grady");
    memory->create(gradz,num_planes*local_x*local_y*local_z,"gradz");
    //theta -- angle between dislocation line tangent and Burges vector
    memory->create(theta,num_planes*local_x*local_y*local_z,"theta");

    memory->create(fcore,num_planes*local_x*local_y*local_z,"fcore");//Check what NP is!
    memory->create(f_core,num_planes*local_x*local_y*local_z,"f_core"); // added for MPI in extended core energy by Claire 7/31/18
    memory->create(df1core,num_planes*local_x*local_y*local_z,"df1core");
    memory->create(df2core,num_planes*local_x*local_y*local_z,"df2core");
    memory->create(df3core,num_planes*local_x*local_y*local_z,"df3core");
    memory->create(dE_core,slip_systems*local_x*local_y*local_z,"dE_core");

    memory->create(xn,slip_systems,dimension,"xn");
    memory->create(xb,slip_systems,dimension,"xb");

    memory->create(data_sigma,2.0*dimension*dimension*local_x*local_y*local_z,"data_sigma");
    memory->create(tau,slip_systems,"tau");
    //memory->create(sigma,dimension,dimension,"sigma");
    //memory->create(deltasig,dimension,dimension,"deltasig");
    memory->create(sigma_rot,dimension,dimension,"sigma_rot");
    memory->create(local_sigma,dimension,dimension,"local_sigma");
    memory->create(ave_sigma,dimension,dimension,"ave_sigma");

    memory->create(eps,slip_systems,dimension,dimension,"eps");
    memory->create(data_eps,2*dimension*dimension*local_x*local_y*local_z,"data_eps");
    memory->create(data_epsd,2*dimension*dimension*local_x*local_y*local_z,"data_epsd");
    memory->create(avepsd,dimension,dimension,"avepsd");
    memory->create(ave_epsd,dimension,dimension,"ave_epsd");
    memory->create(ave_eps,dimension,dimension,"ave_eps");
    memory->create(avepst,dimension,dimension,"avepst");
    memory->create(avepsts,dimension,dimension,"avepsts");
    memory->create(aveps,dimension,dimension,"aveps");

    memory->create(delta,num_planes*local_x*local_y*local_z,"delta");
    memory->create(ddelta,num_planes*local_x*local_y*local_z,"ddelta");

    for(int i=0; i<dimension; i++){
      for(int j=0; j<dimension; j++){
        local_sigma[i][j] = 0;
        ave_sigma[i][j] = 0;
        avepsd[i][j] = 0;
        ave_epsd[i][j] = 0;
        ave_eps[i][j] = 0;
        avepst[i][j] = 0;
        avepsts[i][j] = 0;
        aveps[i][j] = 0;
      }
    }
  }
  /* ----------------------------------------------------------------------
  Calculate sigma_rot
  ------------------------------------------------------------------------- */

  void FFTW_Slab::rotate_stress()
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
          for(k=0; k<ND; k++)
          for(l=0; l<ND; l++){
            sigma_rot[i][j] += prim[i][k]*sigma[k][l]*invprim[l][j];
          }
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

  void FFTW_Slab::frec()
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

void FFTW_Slab::stiffness()
{
  double Caux[dimension][dimension][dimension][dimension];
  double C44 = material->C44;
  double C12 = material->C12;
  double C11 = material->C11;
  /* set Cijkl*/

  for (int i=0; i<dimension; i++) {
    for (int j=0; j<dimension; j++) {
      for (int k=0; k<dimension; k++) {
        for (int m=0; m<dimension; m++) {
          C[i][j][k][m] = C44 * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+C12*DELTA(i,j)*DELTA(k,m);

          /*printf("in BB %d %d %d %d %d\n",dimension, i,j,k,m);*/
        }
      }
    }
  }
  // Rotate Cijkl
  // for (int i=0; i<dimension; i++) {
  //   for (int j=0; j<dimension; j++) {
  //     for (int k=0; k<dimension; k++) {
  // 	for (int m=0; m<dimension; m++) {
  // 	  C[i][j][k][m] = 0;
  // 	  for (int n=0; n<dimension; n++) {
  // 	    for (int u=0; u<dimension; u++) {
  // 	      for (int v=0; v<dimension; v++) {
  // 		for (int w=0; w<dimension; w++) {
  // 		  C[i][j][k][m] += four[i][n]*four[j][u]*four[k][v]*four[m][w]*Caux[n][u][v][w];
  // 		}
  // 	      }
  // 	    }
  // 	  }
  // 	}
  //     }
  //   }
  // }
}

/* ----------------------------------------------------------------------
Calculate Green's function
------------------------------------------------------------------------- */

void FFTW_Slab::greens_function()
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

void FFTW_Slab::Bmatrix()
{
  int i, j, k, l, m, n, u, v, k1, k2, k3, ka, kb, nv, nb, nfreq;
  int is, js, ks;
  double fkr;
  //double C[dimension][dimension][dimension][dimension];
  //double Crot[dimension][dimension][dimension][dimension];
  double A[dimension][dimension][dimension][dimension];

  double B[slip_systems][slip_systems][local_x][local_y][local_z];
  double G[dimension][dimension];
  //double Grot[dimension][dimension];
  double fk[dimension];
  double xn_temp[slip_systems][dimension], xb_temp[slip_systems][dimension], eps_temp[slip_systems][dimension][dimension];
  double xnu, fk2, fk4, fka, fkb, ll, mu, young;
  double C44 = material->C44;
  double C12 = material->C12;
  double C11 = material->C11;
  double dslip = material->dslip;
  // int flag_aniso = materials->flag_aniso;

  mu = C44-(2.0*C44+C12-C11)/5.0;
  ll = C12-(2.0*C44+C12-C11)/5.0;
  young = mu*(3*ll+2*mu)/(ll+mu);
  xnu = young/2.0/mu-1.0;

  material->mu = mu;
  material->ll = ll;
  material->young = young;
  material->xnu = xnu;

  //xnu = nu; //C12/(2.0*(C44+C12));

  if(me == 0){
    printf("Bmatrix mu %lf, ll %lf, young %lf, nu %lf\n", mu, ll, young, xnu);
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
  // /* set Cijkl*/

  // for (i=0; i<dimension; i++) {
  //   for (j=0; j<dimension; j++) {
  //     for (k=0; k<dimension; k++) {
  // 	for (m=0; m<dimension; m++) {
  // 	  C[i][j][k][m] = C44 * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+C12*DELTA(i,j)*DELTA(k,m);
  // 	  Crot[i][j][k][m] = 0.0;
  // 	  A[i][j][k][m] = 0.0;
  // 	  Arot[i][j][k][m] = 0.0;
  // 	  /*printf("in BB %d %d %d %d %d\n",dimension, i,j,k,m);*/
  // 	}
  //     }
  //   }
  // }

  // // Rotate Cijkl
  // for (i=0; i<dimension; i++) {
  //   for (j=0; j<dimension; j++) {
  //     for (k=0; k<dimension; k++) {
  // 	for (m=0; m<dimension; m++) {
  // 	  for (n=0; n<dimension; n++) {
  // 	    for (u=0; u<dimension; u++) {
  // 	      for (v=0; v<dimension; v++) {
  // 		for (w=0; w<dimension; w++) {
  // 		  Crot[i][j][k][m] = four[i][n]*four[j][u]*four[k][v]*four[m][w]*C[n][u][v][w];
  // 		}
  // 	      }
  // 	    }
  // 	  }
  // 	}
  //     }
  //   }
  // }



  for (i=0; i<slip_systems; i++) {
    for (j=0; j<dimension; j++) {
      xn_temp[i][j]=xn[i][j];
      xb_temp[i][j]=xb[i][j];
    }
  }

  /*set eps */
  // DO WE NEED THESE TWO LOOPS???

  for (i=0; i<dimension; i++) {
    for (j=0; j<dimension; j++) {
      for (k=0; k<slip_systems;k++){
        eps_temp[k][i][j] = xb_temp[k][i]*xn_temp[k][j]/dslip;
        eps[k][i][j] = eps_temp[k][i][j];
      }
    }
  }

  // for (i=0; i<dimension; i++) {
  //   for (j=0; j<dimension; j++) {
  //     for (k=0; k<slip_systems;k++){
  // 	eps[k][i][j] = eps_temp[k][i][j];
  //     }
  //   }
  // }

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
        if(fk2>0){
          for (m=0; m<dimension; m++) {
            for (n=0; n<dimension; n++) {
              for (u=0; u<dimension; u++) {
                for (v=0; v<dimension; v++) {
                  A[m][n][u][v] = 0.0;

                  for(i=0; i<dimension; i++) {
                    for(j=0; j<dimension; j++) {
                      for(k=0; k<dimension; k++) {

                        G[k][i] = (2.0 * DELTA(i,k)/fk2-1.0/(1.0-xnu)*fk[i]*fk[k]/fk4)/(2.0*C44);
                        for(l=0; l<dimension; l++) {
                          A[m][n][u][v] = A[m][n][u][v] - C[k][l][u][v]*C[i][j][m][n]*G[k][i]*fk[j]*fk[l] ;
                        }
                      }
                    }
                  }
                  A[m][n][u][v] = A[m][n][u][v]+C[m][n][u][v];
                  /*printf("A %d %d %d %d %f %f %f %f\n",m,n,u,v,fk[0],fk[1],fk[2],A[m][n][u][v]);*/
                }
              }
            }
          }

        } /*if fk2 */
        for(ka=0;ka<slip_systems;ka++){
          for(kb=0;kb<slip_systems;kb++){
            B[ka][kb][k1][k2][k3] = 0.0;
            for (m=0; m<dimension; m++){
              for (n=0; n<dimension; n++) {
                for (u=0; u<dimension; u++) {
                  for (v=0; v<dimension; v++) {
                    B[ka][kb][k1][k2][k3]= B[ka][kb][k1][k2][k3] + A[m][n][u][v]*eps_temp[ka][m][n]*eps_temp[kb][u][v];
                  }
                }
              }
            }

            nb = nfreq +(ka)*local_x*local_y*local_z+(kb)*local_x*local_y*local_z*slip_systems;
            BB[nb] = B[ka][kb][k1][k2][k3]/mu;
            /*printf("%lf %lf %lf %lf \n", fx[nfreq], fy[nfreq], fz[nfreq], BB[nb]);*/
          } /*ka*/
        }/* kb*/


      }	/*k1*/
    }	/*k2*/
  }	/*k3*/


  return;
}

/* ----------------------------------------------------------------------
generate FF matrix
------------------------------------------------------------------------- */

void FFTW_Slab::Fmatrix()
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

  mu = material->mu;
  ll = material->ll;
  young = material->young;
  xnu = material->xnu;

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
  }/*k1,k2,k3*/
  return;
}
/* ----------------------------------------------------------------------
resolve shear stress
------------------------------------------------------------------------- */

void FFTW_Slab::resolSS_Schmid()
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

void FFTW_Slab::resolSS_non_Schmid()
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

void FFTW_Slab::fft_send(void *buf, int count, MPI_Datatype datatype, int dest,
  int tag, MPI_Comm comm)
  {
    MPI_Send(buf,count,datatype,dest,tag,comm);
  }

  /* ----------------------------------------------------------------------
  MPI receive function
  ------------------------------------------------------------------------- */

  void FFTW_Slab::fft_recv(void *buf, int count, MPI_Datatype datatype, int dest,
    int tag, MPI_Comm comm, MPI_Status *status)
    {
      MPI_Recv(buf,count,datatype,dest,tag,comm,status);
    }
    /* ----------------------------------------------------------------------
    grad(xi) and line character angle theta
    ------------------------------------------------------------------------- */
    void FFTW_Slab::gradient()
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
      double xi_r[nmpi], xi_l[nmpi], mpi_r[nmpi], mpi_l[nmpi]; //order parameters in the skin region
      double d1,d2,d3;
      d1 = d2 = d3 = 1.0;
      double tang[3], mag;
      MPI_Status status;

      //initialize
      for (i=0; i<nmpi; i++) {
        xi_r[i] = 0.0;
        xi_l[i] = 0.0;
        mpi_r[i] = 0.0;
        mpi_l[i] = 0.0;
      }
      for (i=0; i<nmpi*lN1; i++) {
        gradx[i] = 0.0;
        grady[i] = 0.0;
        gradz[i] = 0.0;
        theta[i] = 0.0;
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

      // transformation of the gradient from the local to global grid
      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          gradx[index] = gradx[index]*four[0][0] + grady[index]*four[1][0] + gradz[index]*four[2][0];
          grady[index] = gradx[index]*four[0][1] + grady[index]*four[1][1] + gradz[index]*four[2][1];
          gradz[index] = gradx[index]*four[0][2] + grady[index]*four[1][2] + gradz[index]*four[2][2];
        }
      }

      double dir;
      for(isa=0;isa<NS;isa++){

        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          tang[0] = grady[index]*xn[isa][2] - gradz[index]*xn[isa][1];
          tang[1] = gradz[index]*xn[isa][0] - gradx[index]*xn[isa][2];
          tang[2] = gradx[index]*xn[isa][1] - grady[index]*xn[isa][0];
          mag = sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]);
          if (mag<=0.1) {
            theta[index]= -1; //define
          }
          else {
            tang[0] /= mag;
            tang[1] /= mag;
            tang[2] /= mag;
            theta[index] = acos(tang[0]*xb[isa][0]+tang[1]*xb[isa][1]+tang[2]*xb[isa][2]);
            dir = (xb[isa][1]*tang[2]-xb[isa][2]*tang[1])*xn[isa][0]+
            (xb[isa][2]*tang[0]-xb[isa][0]*tang[2])*xn[isa][1]+
            (xb[isa][0]*tang[1]-xb[isa][1]*tang[0])*xn[isa][2];
            if (dir<0.0) {
              theta[index]=acos(-1)-theta[index];
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
    }

    /* -----------------------------------------------------------------------
    bcc core energy for perfect dislocations
    ---------------------------------------------------------------------*/
    void FFTW_Slab::core_energy_bcc_perfect()
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
    void FFTW_Slab::core_energy_perfect()
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
      // Checking without core energy
      //An = 0;

      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
        for(k=0;k<N3;k++){
          index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
          E_core += An*(sin(M_PI*data_fftw[index].re)*sin(M_PI*data_fftw[index].re))/N1;   //To make it general fft->data_fftw hasz to be general
          dE_core[index] = An*M_PI*sin(2.0*M_PI*data_fftw[index].re);
        }/*ijk*/
      }/*isa*/
    }

    /* -----------------------------------------------------------------------
    core energy for variable USFE
    ---------------------------------------------------------------------*/
    void FFTW_Slab::core_energy_mpea()
    {
      int i, j, k, isa, index, num, core_index;
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

      double * A_mpea = material->A_mpea;
      double b = material->b;

      double mu = material->mu;
      double ll = material->ll;
      double young = material->young;
      double xnu = material->xnu;

      E_core = 0;

      for(isa=0;isa<NS;isa++){
        for(i=0;i<lN1;i++)
          for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
              index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
              core_index = (lxs+i)*N2*N3 + j*N3 + k + isa*N1*N2*N3;
              E_core += A_mpea[core_index]/(mu*dslip*b)*(sin(M_PI*data_fftw[index].re)*sin(M_PI*data_fftw[index].re))/N1;
              dE_core[index] = A_mpea[core_index]/(mu*dslip*b)*M_PI*sin(2.0*M_PI*data_fftw[index].re);
        }/*ijk*/
      }/*isa*/
    }

    /* -----------------------------------------------------------------------
    core energy extended dislocations
    ---------------------------------------------------------------------*/
    /* Claire added extended energy 07/30/18 to account for partial dislocations
    and the Schoeck parameterization. Can be used for both fcc {111} planes and
    hcp {0001} basal planes to model partial dislocations*/
    void FFTW_Slab::core_energy_extended()
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
    core energy sine approximation
    ---------------------------------------------------------------------*/
    /* Added by Claire 8/1/18: a 1D Sine approximations for use modeling hcp
    non-basal slip planes with GSFE that suggests partial dislocations */
    void FFTW_Slab::core_energy_sine()
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
    void FFTW_Slab::core_energy_pyrII()
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

    void FFTW_Slab::core_energy_USFE_angle_tau()
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
    void FFTW_Slab::project_core_energy()
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
        Initialize inner loop
        ---------------------------------------------------------------------*/
        void FFTW_Slab::init_loop()
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
        void FFTW_Slab::internal_energy()
        {
          int lN1 = local_x;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;
          int index=0, index2=0, nb=0;

          for (int i=0; i<lN1*N2*N3*NS; i++){
            data_real[i].re = 0;
            data_real[i].im = 0;
          }
          for(int isa=0;isa<NS;isa++){
            for(int isb=0;isb<NS;isb++){
              for(int i=0;i<lN1;i++)
              for(int j=0;j<N2;j++)
              for(int k=0;k<N3;k++){
                index  = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
                index2 = i*N2*N3 + j*N3 + k + isb*lN1*N2*N3;
                nb     = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3 + isb*lN1*N2*N3*NS;
                data_real[index].re += data_fftw[index2].re * BB[nb];
                data_real[index].im += data_fftw[index2].im * BB[nb];
              }
            }
          }
        }

        /* -----------------------------------------------------------------------
        Update order parameter
        ---------------------------------------------------------------------*/
        void FFTW_Slab::update_order_parameter()
        {
          int lN1 = local_x;
          int N1 = nx;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;
          int index=0, na=0, na0=0, na1=0;
          double xinormlocal=0.0;           // stores the norm of the increment of xi
          double xirep=0.0, xiimp=0.0;      // store previous re and im values of xi
          double xi_ave=0.0;                // Local average of the order parameter
          int nsize = N1*N2*N3;

          for(int isa=0;isa<NS;isa++){
            for(int i=0;i<lN1;i++){
              for(int j=0;j<N2;j++){
                for(int k=0;k<N3;k++){
                  na0 = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
                  index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
                  na = 2*(i*N2*N3 + j*N3 + k);
                  na1 = na0+1;

                  //Ginzburg-Landau Equation for real and imag parts

                  if(xo[index] == 0.0){
                    xirep = xi[0][na0];
                    xiimp = xi[0][na1];
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
          MPI_Allreduce(&xi_ave, &xiave, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(&xinormlocal, &xinorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          xinorm = sqrt(xinorm);
        }
        /* ------------------------------------------------------------------
        Prepares the next iteration
        ------------------------------------------------------------------ */
        void FFTW_Slab::prepare_next_itr()
        {
          int lN1 = local_x;
          int N1 = nx;
          int N2 = local_y;
          int N3 = local_z;
          int NS = slip_systems;
          int index=0, na=0, na0=0, na1=0;

          for(int isa=0;isa<NS;isa++){
            for(int i=0;i<lN1;i++){
              for(int j=0;j<N2;j++){
                for(int k=0;k<N3;k++){
                  na0 = 2*(i*N2*N3 + j*N3 + k + isa*lN1*N2*N3);
                  index = i*N2*N3 + j*N3 + k + isa*lN1*N2*N3;
                  na = 2*(i*N2*N3 + j*N3 + k);
                  na1 = na0+1;
                  data_fftw[index].re = xi[0][na0];
                  data_fftw[index].im = xi[0][na1];

                }
              }
            }
          }
        }
        /* -----------------------------------------------------------------------
        temporary stores data_fftw
        ---------------------------------------------------------------------*/
        void FFTW_Slab::temp_stor_data()
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
        void FFTW_Slab::reset_data()
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
        void FFTW_Slab::prep_forward()
        {
          if(mode == 1)
          forward_mode1();
          else if(mode == 2)
          forward_mode2();
        }

        /* ------------------------------------------------------------------
        Prepares the FFT backward
        ------------------------------------------------------------------ */
        void FFTW_Slab::prep_backward()
        {
          if(mode == 1)
          backward_mode1();
          else if(mode == 2)
          backward_mode2();
        }


        /* -----------------------------------------------------------------------
        forward FFT
        ---------------------------------------------------------------------*/
        void FFTW_Slab::forward_mode1()
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
        void FFTW_Slab::forward_mode2()
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
        backward FFT
        ---------------------------------------------------------------------*/
        void FFTW_Slab::backward_mode1()
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
        void FFTW_Slab::backward_mode2()
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
        Initial Crystal to nothing
        ---------------------------------------------------------------------*/
        void FFTW_Slab::initial_sxtal_NoLoop()
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
        void FFTW_Slab::initial_sxtal_mode1()
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
        void FFTW_Slab::initial_sxtal_mode2()
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

        void FFTW_Slab::initial_sxtal_1LOrthoFCC()
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
        void FFTW_Slab::initial_sxtal_1LNonOrthoFCC()
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
        void FFTW_Slab::initial_sxtal_1LInclinedFCC()
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
        void FFTW_Slab::initial_sxtal_1LOrthoBCC()
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
        void FFTW_Slab::initial_sxtal_2LOrtho()
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
        void FFTW_Slab::initial_sxtal_StraightOrtho()
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
        	  xi[0][na1] = 0.0;
        	  xo[index] = 0.0;

            if(((k-local_z/2.0)<eps) &&     // {001} plane
              ((k-local_z/2.0)>-eps) &&
              (j>(N2/4.0)) && (j<(3.0*N2/4.0)) && //((lxs+i)>(N1/4.0)) && ((lxs+i)<(3.0*N1/4.0)) && //
              (is==0)){


              xi[0][na0]=1.0;
              xi_sum[0][na] = xi_sum[0][na] + xi[0][na0];

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
        void FFTW_Slab::initial_sxtal_StraightNonOrtho()
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
        void FFTW_Slab::initial_sxtal_3D3ShcpBasal()
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
        void FFTW_Slab::initial_sxtal_3D3ShcpNotch()
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

      void FFTW_Slab::initial_sxtal_frank_read()
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

  void FFTW_Slab::initial_sxtal_3SBCC()
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
      data_fftw[index].re = 0.0;
      data_fftw[index].im = 0.0;

      if(is == 0 && (lxs+i) == 0 && k >= N3/4 && k < N3*3/4){
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

  void FFTW_Slab::initial_sxtal_loop()
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

  void FFTW_Slab::initial_sxtal_edge()
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
  void FFTW_Slab::stressfree_strain()
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
    int i, j, k, l, is, nb, ida, idb;
    double S[ND][ND][ND][ND];
    double S44 = material->S44;
    double S12 = material->S12;
    double mu = material->mu;

    /*set S matrix*/
    for (i=0; i<ND; i++)
    for (j=0; j<ND; j++)
    for (k=0; k<ND; k++)
    for (l=0; l<ND; l++)
    {
      S[i][j][k][l] = S44/4 * (DELTA(i,k)*DELTA(j,l)+DELTA(i,l)*DELTA(j,k))+S12*DELTA(i,j)*DELTA(k,l);
    }


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
void FFTW_Slab::total_average_strain()
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
    }
  }
  for(ida=0;ida<ND;ida++){
    for (idb=0;idb<ND;idb++){
      for(i=0;i<lN1;i++){
        for(j=0;j<N2;j++){
          for (k=0;k<N3;k++){
            na0 = 2*(k + j*N3 + i*N2*N3 + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
            aveps[ida][idb] += data_eps[na0]; //aveps only appears here to get total average strain in the box
          }
        }
      }
      MPI_Reduce(&aveps[ida][idb], &ave_eps[ida][idb], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      ave_eps[ida][idb] /= nsize;
    } //idb
  } //ida

}

/* -----------------------------------------------------------------------
strain
---------------------------------------------------------------------*/
void FFTW_Slab::strain()
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
  int i, j, k, l, is, na0, na1, nb, psys, ida, idb, index, index2;
  int na11, na12, na13, na21, na22, na23, na31, na32, na33, ia, ib;
  fftwnd_mpi_plan iiplan;

  iiplan = fftw3d_mpi_create_plan(world, N1, N2, N3, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (i=0; i<lN1*N2*N3*ND*ND; i++){
    data_strain[i].re = 0.0;
    data_strain[i].im = 0.0;
  }
  for (i=0; i<lN1*N2*N3; i++){
    work_strain[i].re = 0.0;
    work_strain[i].im = 0.0;
  }
  for (i=0; i<2*lN1*N2*N3*ND*ND; i++){
    data_eps[i] = 0.0;
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

  for(i=0;i<ND;i++){
    for (j=0;j<ND;j++){
      psys = i*lN1*N2*N3 + j*lN1*N2*N3*ND;
      fftwnd_mpi(iiplan, 1, data_strain+psys, work_strain, FFTW_NORMAL_ORDER); /* Inverse FFT (multiple)*/
    }
  }
  for (i=0; i<lN1*N2*N3*ND*ND; i++){
    data_strain[i].re = data_strain[i].re/(nsize);
    data_strain[i].im = data_strain[i].im/(nsize);
  }

  // calculate stress free strain
  stressfree_strain();
  //add in other two terms in strain (added already into avepst: epsilon bar 0 and S x sigma)
  for(ida=0;ida<ND;ida++){
    for (idb=0;idb<ND;idb++){
      for(i=0;i<lN1;i++){
        for(j=0;j<N2;j++){
          for (k=0;k<N3;k++){
            index = k + j*N3 + i*N2*N3 + ida*lN1*N2*N3 + idb*lN1*N2*N3*ND;
            data_strain[index].re += avepst[ida][idb];     // Adding epsilon bar Wang 2002
          }
        }
      }
    }
  }

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
void FFTW_Slab::stress()
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
  double mu, xnu, young, ll;

  double C44 = material->C44;
  double C12 = material->C12;
  double C11 = material->C11;

  // Why do we have to redefine this?????

  mu = C44-(2.0*C44+C12-C11)/5.0;
  ll = C12-(2.0*C44+C12-C11)/5.0;
  young = mu*(3*ll+2*mu)/(ll+mu);
  xnu = young/2.0/mu-1.0;

  for (i=0; i<2*lN1*N2*N3*ND*ND; i++){
    data_sigma[i] =0;
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
                data_sigma[na] += C[ida][idb][m][l]*(data_strain[index].re - data_epsd[na0] - avepsts[m][l]);
                data_sigma[na+1] = 0.0;
              }
            }
            //data_sigma[na]+=sigma[ida][idb]*mu;
            data_sigma[na]+=sigma_rot[ida][idb]*mu;
            local_sigma[ida][idb] += data_sigma[na];
          }
        }
      }
      MPI_Allreduce(&local_sigma[ida][idb], &ave_sigma[ida][idb], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      ave_sigma[ida][idb] /= nsize;
    }
  }
}
