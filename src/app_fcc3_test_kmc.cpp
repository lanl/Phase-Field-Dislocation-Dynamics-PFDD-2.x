/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
   ------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_fcc3_test_kmc.h"
#include "fft.h"
#include "material.h" //added this 
#include "solve.h"
#include "random_park.h"
#include "error.h"
#include "memory.h"

using namespace PFDD_NS;;

/* ---------------------------------------------------------------------- */

AppFCC3TestKMC::AppFCC3TestKMC(PFDD_C *pfdd_p, int narg, char **arg) :
  App(pfdd_p,narg,arg)
{
  ndouble = 0;
  timestep = 1.0;
  slip_systems = 3;  // The # of slip systems has to be allocated in the constructure
                      //  fot the FFT class to be able to allocate arrays
  num_planes = 1;  // The # of plane normals has to be allocated in the constructure
                      //  fot the FFT class to be able to allocate arrays

  CD = 0.25;           // Dislocation mobility
  dflag = 1;

  NE = 1;             // total # of events for KMC
  // memory->create(event_cond,NE,"event_cond");
  memory->create(event_time,NE,"event_time");
  memory->create(pre_fac,NE,"pre_fac");
  memory->create(E_act,NE,"E_act");
  memory->create(rate,NE,"rate");
  memory->create(stress_dependence,NE,"stress_dependence"); // NEW LINE 
  memory->create(vol_ac,NE,dimension,dimension,"vol_ac");
/*
  for(int i=0; i<NE; i++){
    event_time[i] = -1.0;
  }
  */
  // cur_time = 0.0;
  // dt_kmc = 0.0;
  beta = 0.0;
  kflag = 0; // default flag for KMC

  if(setup_flag == -1)
    setup_flag = 3; //default to a screw dipole
  if(core_flag == -1)
    core_flag = 1; //default to sin^2
}

/* ---------------------------------------------------------------------- */

AppFCC3TestKMC::~AppFCC3TestKMC()
{
  delete [] sites;
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
   ------------------------------------------------------------------------- */

void AppFCC3TestKMC::grow_app()
{
  darray = fft->xi;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
   ---------------------------------------------------------------------- */

void AppFCC3TestKMC::init_app()
{
  delete [] sites;
  sites = new int[1 + maxneigh];

  grow_app();
  set_slip();
  initial_sxtal();

  stopsteps = static_cast<int>(floor(stoptime/timestep));

  // char temp1[500], temp2[500];

  // sprintf(temp1,"dumps/slip_%d.txt",me);
  // sprintf(temp2,"dumps/energy_%d.txt",me);
  // sfile = fopen(temp1,"w");
  // efile = fopen(temp2, "w");
  // fccflag = 0;
}
/* ----------------------------------------------------------------------
   input specific for this app
   ------------------------------------------------------------------------- */
void AppFCC3TestKMC::input_app(char *command, int narg, char **arg)
{

  if (strcmp(command,"events") == 0) {
    if (narg < 10) error->all(FLERR,"Illegal app_style command");
    int iarg=0;
    int ievent  = atoi(arg[iarg]);
    iarg++;
    pre_fac[ievent] = atoi(arg[iarg]);
    iarg++;
    E_act[ievent] = atoi(arg[iarg]);
    iarg++;
    vol_ac[ievent][0][0] = atoi(arg[iarg]);
    iarg++;
    vol_ac[ievent][1][1] = atoi(arg[iarg]);
    iarg++;
    vol_ac[ievent][2][2] = atoi(arg[iarg]);
    iarg++;
    vol_ac[ievent][0][1] = vol_ac[ievent][1][0] = atoi(arg[iarg]);
    iarg++;
    vol_ac[ievent][1][2] = vol_ac[ievent][2][1] = atoi(arg[iarg]);
    iarg++;
    vol_ac[ievent][2][0] = vol_ac[ievent][0][2] = atoi(arg[iarg]);
    iarg++; //NEW LINE
    stress_dependence[ievent] = atoi(arg[iarg]);
    if(me == 0){
      if (logfile) fprintf(logfile,"event No %d, stress_dependence = %lf\n", ievent, stress_dependence[ievent]);
      if (screen) fprintf(screen,"event No %d, stress_dependence = %lf\n", ievent, stress_dependence[ievent]);
    }
    // if(me == 0){
    //   if (logfile) fprintf(logfile,"event No %d, pre_fac = %lf, E_act = %lf\n", ievent, pre_fac[ievent], E_act[ievent]);
    //   if (screen) fprintf(screen,"event No %d, pre_fac = %lf, E_act = %lf\n", ievent, pre_fac[ievent], E_act[ievent]);
    // }
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
   ------------------------------------------------------------------------- */

// double AppFCC::site_energy(int i)
// {
//   int isite = spin[i];
//   int eng = 0;
//   for (int j = 0; j < numneigh[i]; j++)
//     if (isite != spin[neighbor[i][j]]) eng++;
//   return (double) eng;
// }

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   flip randomly to either spin
   null bin extends to size 2
   ------------------------------------------------------------------------- */

// void AppFCC::site_event_rejection(int i, RandomPark *random)
// {
//   int oldstate = spin[i];
//   double einitial = site_energy(i);

//   // event = random spin from 1 to 2, including self

//   if (random->uniform() < 0.5) spin[i] = 1;
//   else spin[i] = 2;
//   double efinal = site_energy(i);

//   // accept or reject via Boltzmann criterion

//   if (efinal <= einitial) {
//   } else if (temperature == 0.0) {
//     spin[i] = oldstate;
//   } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
//     spin[i] = oldstate;
//   }

//   if (spin[i] != oldstate) naccept++;

//   // set mask if site could not have changed
//   // if site changed, unset mask of sites with affected propensity
//   // OK to change mask of ghost sites since never used

//   if (Lmask) {
//     if (einitial < 0.5*numneigh[i]) mask[i] = 1;
//     if (spin[i] != oldstate)
//       for (int j = 0; j < numneigh[i]; j++)
// 	mask[neighbor[i][j]] = 0;
//   }
// }

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
   ------------------------------------------------------------------------- */

// double AppFCC::site_propensity(int i)
// {
//   // event = spin flip

//   int oldstate = spin[i];
//   int newstate = 1;
//   if (oldstate == 1) newstate = 2;

//   // compute energy difference between initial and final state
//   // if downhill or no energy change, propensity = 1
//   // if uphill energy change, propensity = Boltzmann factor

//   double einitial = site_energy(i);
//   spin[i] = newstate;
//   double efinal = site_energy(i);
//   spin[i] = oldstate;

//   if (efinal <= einitial) return 1.0;
//   else if (temperature == 0.0) return 0.0;
//   else return exp((einitial-efinal)*t_inverse);
// }

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
   ------------------------------------------------------------------------- */

// void AppFCC::site_event(int i, RandomPark *random)
// {
//   int m;

//   // perform event = spin flip

//   if (spin[i] == 1) spin[i] = 2;
//   else spin[i] = 1;

//   // compute propensity changes for self and neighbor sites
//   // ignore update of neighbor sites with isite < 0

//   int nsites = 0;
//   int isite = i2site[i];
//   sites[nsites++] = isite;
//   propensity[isite] = site_propensity(i);

//   for (int j = 0; j < numneigh[i]; j++) {
//     m = neighbor[i][j];
//     isite = i2site[m];
//     if (isite < 0) continue;
//     sites[nsites++] = isite;
//     propensity[isite] = site_propensity(m);
//   }

//   solve->update(nsites,sites,propensity);
// }

/* ----------------------------------------------------------------------
   New
   returns the spin of one site
   ------------------------------------------------------------------------- */

// int AppFCC::site_spin (int i) {
//   if ((i<0)||(i>nlocal))
//     error->all (FLERR, "Invalid site index");
//   return spin[i];
// }


void AppFCC3TestKMC::set_slip()
{
  int i, j, k;
  double s3;
  double s2;
  double s6;
  double xnf=0,ynf=0,znf=0;
  double xbf=0,ybf=0,zbf=0;

  s3 = sqrt(3);
  s2 = sqrt(2);
  for (i=0;i<3 ; i++) {
    fft->xn[i][0]= 1.0/s3;
    fft->xn[i][1]= 1.0/s3;
    fft->xn[i][2]= 1.0/s3;
  }
  // for (i=3;i<6 ; i++) {
  //   fft->xn[i][0]= -1.0/s3;
  //   fft->xn[i][1]= -1.0/s3;
  //   fft->xn[i][2]= 1.0/s3;
  // }
// MUST BE IN THIS ORDER FOR LOOKUP TABLE: <110> <101> <011>

  fft->xb[0][0] = -1.0/s2;
  fft->xb[0][1] = 1.0/s2;
  fft->xb[0][2] = 0.0;

  fft->xb[1][0] = 1.0/s2;
  fft->xb[1][1] = 0.0;
  fft->xb[1][2] = -1.0/s2;

  fft->xb[2][0] = 0.0;
  fft->xb[2][1] = -1.0/s2;
  fft->xb[2][2] = 1.0/s2;

  // fft->xb[3][0] = -1.0/s2;
  // fft->xb[3][1] = 1.0/s2;
  // fft->xb[3][2] = 0.0;

  // fft->xb[4][0] = 1.0/s2;
  // fft->xb[4][1] = 0.0;
  // fft->xb[4][2] = 1.0/s2;

  // fft->xb[5][0] = 0.0;
  // fft->xb[5][1] = -1.0/s2;
  // fft->xb[5][2] = -1.0/s2;

}
/* ----------------------------------------------------------------------
   Check if KMC events happen and their time
   ------------------------------------------------------------------------- */

void AppFCC3TestKMC::events()
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N1 = fft->nx;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int index, i, j, k;
  /*
  int event_sum = 0;
  for(int i=0; i<NE; i++){
    event_cond[i] = event_check(i);
    event_sum += event_cond[i];
  }
  if (event_sum == 0) return;
  for(int i=0; i<NE; i++){
    if (event_cond[i] == 1 && event_time[i] < 0.0){
      event_time[i] == calculate_event_time(i);
    }
  }
  */
  int is = 0; // slip system 1
  int na0, na1, index0, index1, index2, index3;
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1 = i*N2*N3 + j*N3 + k + (3*is+0)*lN1*N2*N3;
    index2 = i*N2*N3 + j*N3 + k + (3*is+1)*lN1*N2*N3;
    index3 = i*N2*N3 + j*N3 + k + (3*is+2)*lN1*N2*N3;
    na0 = 2*(i*N2*N3+j*N3+k+is*lN1*N2*N3);
    na1 = na0+1;
    if (((lxs+i) == 0 && ((k == int(N3/8.0)) || (k == int(N3*7.0/8.0))))){
    if (fft->picked_d[index1]>0.0 
       //|| (j >= 49.9 && j<= 50.1)||(j >= 14.9 && j<= 15.1))
       ) {
     /* if(me == 0){
        if (logfile) fprintf(logfile,"Meet %d, %d, %d\n", i,j,k);
        if (screen) fprintf(screen,"Meet %d, %d, %d\n", i,j,k);
      }*/
      if (fft->xo[index1] <= 0.0) {
        fft->xo[index1] = 1.0;
        fft->xo[index2] = 1.0;
        fft->xo[index3] = 1.0;
     /*   if(me == 0){
          if (logfile) fprintf(logfile,"xo %lf\n", fft->xo[index]);
          if (screen) fprintf(screen,"xo %lf\n", fft->xo[index]);
        }*/
        calculate_event_time(index1,j);
        // fft->picked_d[index1] = dt_kmc;
        barrier(j,k,1);
      }
      else if (fft->xo[index1] == 1.0){
        if (fft->picked_d[index1]<=1.0) {
         // fft->xo[index] = -15.0;
          fft->xo[index1] = -1.0;// * fft->xi[0][na0];
          fft->xo[index2] = -1.0;// * fft->xi[0][na0];
          fft->xo[index3] = -1.0;// * fft->xi[0][na0];
          fft->picked_d[index1] = 1.0;

          if (stress_dependence[0] == 1.0) {
            kmc_stress_update(index1,j);
            if(me == 0){
              if (logfile) fprintf(logfile,"stress dependence activated");
              if (screen) fprintf(screen,"stress dependence activated");
            }
          }
          
          barrier(j,k,0);
        }
        else{
          fft->picked_d[index1] -= 1.0;
        }
      }
    }
      else {
      fft->xo[index1]==0.0;
        }
      // fft->events(index);
    }
  }

}

void AppFCC3TestKMC::kmc_stress_update(int index, int j) // local KMC re-computation 
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N1 = fft->nx;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int ND = dimension;
  double * Hv = fft->material->Hv; // Pull these from the read_activation_file
  double * Pv = fft->material->Pv;

  double kappa = 8.6173E-5; // Boltzmann constant
  beta = temperature * kappa;

  int index1, index2, j1, j2; 
  // double dE = 0.4; // decrement of activation energy 
  double minE = 0.0; // minimum activation energy

  // Activation volume tensor components defined here, calibrated beforehand in OVITO 
  // double vxx = 0.0; //1.0E-5;
  // double vyy = 0.0; //1.0E-5;
  // double vzz = 0.0; //1.0E-5; 
  // double vxy = 0.0; 
  // double vxz = 0.0; 
  // double vyz = 0.0; 
  double vxx = vol_ac[0][0][0];
  double vyy = vol_ac[0][1][1];
  double vzz = vol_ac[0][2][2]; 
  double vxy = vol_ac[0][0][1];
  double vxz = vol_ac[0][1][2];
  double vyz = vol_ac[0][2][0]; 

  double unitconv = 6.2415e-12; //Pa --> eV/ang^3; //Pa --> eV/ang^3
                            
  // indices to left and right of activated grid point in 2-D (flattened) / 1-D
  index1 = index - N3;
  index2 = index + N3;
  j1 = j-1;
  j2 = j+1;

  // PBCs for flattened plane coordinates
  if (index1 < 0) {
    index1 = index + (N2-1)*N3;
  } 
  if (index2 >= lN1*N2*N3) {
    index2 = index - (N2-1)*N3;
  } 
  // PBCs for local interface coordinates
  if (j1 < 0) {
    j1 = j + (N2-1);
  }
  if (j2 >= N2) {
    j2 = j - (N2-1);
  }

  if(me == 0){
    if (logfile) fprintf(logfile,"vxx %.7f, vyy = %.7f, vzz = %.7f,fft->picked[index1]=%.7f\n", vxx, vyy, vzz, fft->picked[index1]);
    if (screen) fprintf(screen,"vxx %.7f, vyy = %.7f, vzz = %.7f,fft->picked[index1]=%.7f\n", vxx, vyy, vzz, fft->picked[index1]);
  }

  if (fft->picked_d[index1]>1.0) { //if grid point has not yet transmitted, then reset KCM with new activation energy
    
    // double pxx1 = compute_stress(index1,1); //define local stresses
    // double pyy1 = compute_stress(index1,2);
    // double pzz1 = compute_stress(index1,3);
    // double pxy1 = compute_stress(index1,4);
    // double pxz1 = compute_stress(index1,5);
    // double pyz1 = compute_stress(index1,6);
    double pxx1 = fft->data_sigma[2*(index1 + 0*lN1*N2*N3 + 0*lN1*N2*N3*ND)]; //define local stresses
    double pyy1 = fft->data_sigma[2*(index1 + 1*lN1*N2*N3 + 1*lN1*N2*N3*ND)];
    double pzz1 = fft->data_sigma[2*(index1 + 2*lN1*N2*N3 + 2*lN1*N2*N3*ND)];
    double pxy1 = fft->data_sigma[2*(index1 + 0*lN1*N2*N3 + 1*lN1*N2*N3*ND)];
    double pxz1 = fft->data_sigma[2*(index1 + 0*lN1*N2*N3 + 2*lN1*N2*N3*ND)];
    double pyz1 = fft->data_sigma[2*(index1 + 1*lN1*N2*N3 + 2*lN1*N2*N3*ND)];

    if(me == 0){
      if (logfile) fprintf(logfile,"local stress at index1:");
      if (logfile) fprintf(logfile,"pxx1 %e, pyy1 %e, pzz1 %e\n",pxx1,pyy1,pzz1);
      if (logfile) fprintf(logfile,"pxy1 %e, pxz1 %e, pyz1 %e\n",pxy1,pxz1,pyz1);
      if (screen) fprintf(screen,"local stress at index1:");
      if (screen) fprintf(screen,"pxx1 %e, pyy1 %e, pzz1 %e\n",pxx1,pyy1,pzz1);
      if (screen) fprintf(screen,"pxy1 %e, pxz1 %e, pyz1 %e\n",pxx1,pyy1,pzz1);
    }

    double dE1 = (pxx1*vxx+pyy1*vyy+pzz1*vzz+2*pxy1*vxy+2*pxz1*vxz+2*pyz1*vyz)*unitconv;

    Hv[j1] = Hv[j1]-dE1;

    if (Hv[j1]<minE){ // check that activation energy does not go below zero or minE
      Hv[j1] = minE;
    }

    for(int i=0;i<NE;i++){
      rate[i] = Pv[j1]*exp(-abs(Hv[j1])/beta);
    } 

    double rannum = ranapp->uniform();
    // dt_kmc = 1 / rate[0] / 2.48 / 1.0E-14 * log(1/rannum);
    // fft->picked_d[index1] = dt_kmc;

    if(me == 0){
      if (logfile) fprintf(logfile,"Hv[j1] %.7f, dE1 %.7f\n", Hv[j1],dE1);
      if (logfile) fprintf(logfile,"resetting dt_kmc to fft->picked[index1]=%.7f\n",fft->picked[index1]);
      if (screen) fprintf(screen,"Hv[j1] %.7f, dE1 %.7f\n", Hv[j1],dE1);
      if (screen) fprintf(screen,"resetting dt_kmc to fft->picked[index1]=%.7f\n",fft->picked[index1]);
    }
  }

  if (fft->picked_d[index2]>1.0) { //if grid point has not yet transmitted, then reset KMC with lower activation energy
    
    // double pxx2 = compute_stress(index2,1); //define local stresses
    // double pyy2 = compute_stress(index2,2);
    // double pzz2 = compute_stress(index2,3);
    // double pxy2 = compute_stress(index2,4);
    // double pxz2 = compute_stress(index2,5);
    // double pyz2 = compute_stress(index2,6);
    double pxx2 = fft->data_sigma[2*(index2 + 0*lN1*N2*N3 + 0*lN1*N2*N3*ND)]; //define local stresses
    double pyy2 = fft->data_sigma[2*(index2 + 1*lN1*N2*N3 + 1*lN1*N2*N3*ND)];
    double pzz2 = fft->data_sigma[2*(index2 + 2*lN1*N2*N3 + 2*lN1*N2*N3*ND)];
    double pxy2 = fft->data_sigma[2*(index2 + 0*lN1*N2*N3 + 1*lN1*N2*N3*ND)];
    double pxz2 = fft->data_sigma[2*(index2 + 0*lN1*N2*N3 + 2*lN1*N2*N3*ND)];
    double pyz2 = fft->data_sigma[2*(index2 + 1*lN1*N2*N3 + 2*lN1*N2*N3*ND)];

    if(me == 0){
      if (logfile) fprintf(logfile,"local stress at index2:");
      if (logfile) fprintf(logfile,"pxx2 %e, pyy2 %e, pzz2 %e\n",pxx2,pyy2,pzz2);
      if (logfile) fprintf(logfile,"pxy2 %e, pxz2 %e, pyz2 %e\n",pxy2,pxz2,pyz2);
      if (screen) fprintf(screen,"local stress at index2:");
      if (screen) fprintf(screen,"pxx2 %e, pyy2 %e, pzz2 %e\n",pxx2,pyy2,pzz2);
      if (screen) fprintf(screen,"pxy2 %e, pxz2 %e, pyz2 %e\n",pxx2,pyy2,pzz2);
    }

    double dE2 = (pxx2*vxx+pyy2*vyy+pzz2*vzz+2*pxy2*vxy+2*pxz2*vxz+2*pyz2*vyz)*unitconv;
    Hv[j2] = Hv[j2]-dE2;

    if (Hv[j2]<minE){ // check that activation energy does not go below zero or minE
      Hv[j2] = minE;
    }

    for(int i=0;i<NE;i++){
      rate[i] = Pv[j2]*exp(-abs(Hv[j2])/beta);
    } 

    double rannum2 = ranapp->uniform();
    // dt_kmc = 1 / rate[0] / 2.48 / 1.0E-14 * log(1/rannum2);
    // fft->picked_d[index2] = dt_kmc;

    if(me == 0){
      if (logfile) fprintf(logfile,"Hv[j2] %.7f, dE2 %.7f\n", Hv[j2],dE2);
      if (logfile) fprintf(logfile,"resetting dt_kmc to fft->picked[index2]=%.7f\n",fft->picked[index2]);
      if (screen) fprintf(screen,"Hv[j2] %.7f, dE2 %.7f\n", Hv[j2],dE2);
      if (screen) fprintf(screen,"resetting dt_kmc to fft->picked[index2]=%.7f\n",fft->picked[index2]);
    }
  }
}

void AppFCC3TestKMC::barrier(int y, int z, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N1 = fft->nx;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int index, i, j, k;
  /*
  int event_sum = 0;
  for(int i=0; i<NE; i++){
    event_cond[i] = event_check(i);
    event_sum += event_cond[i];
  }
  if (event_sum == 0) return;
  for(int i=0; i<NE; i++){
    if (event_cond[i] == 1 && event_time[i] < 0.0){
      event_time[i] == calculate_event_time(i);
    }
  }
  */
  int is = 0; // slip system 1
  int na0, index0, index1, index2, index3;
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index1 = i*N2*N3 + j*N3 + k + (3*is+0)*lN1*N2*N3;
    index2 = i*N2*N3 + j*N3 + k + (3*is+1)*lN1*N2*N3;
    index3 = i*N2*N3 + j*N3 + k + (3*is+2)*lN1*N2*N3;
    if (z < N3/2){
          if ((lxs+i) == 0 && j == y && (k < int(N3/8.0)))
       //|| (j >= 49.9 && j<= 50.1)||(j >= 14.9 && j<= 15.1))
       {
     /* if(me == 0){
        if (logfile) fprintf(logfile,"Meet %d, %d, %d\n", i,j,k);
        if (screen) fprintf(screen,"Meet %d, %d, %d\n", i,j,k);
      }*/
      if (flag == 1  /*&& fft->firstpass[j][0]==0 */) {
        fft->xo[index1] = 1.0;
        fft->xo[index2] = 1.0;
        fft->xo[index3] = 1.0;
     /*   if(me == 0){
          if (logfile) fprintf(logfile,"xo %lf\n", fft->xo[index]);
          if (screen) fprintf(screen,"xo %lf\n", fft->xo[index]);
        }*/
        // fft->picked[index1] = dt_kmc;
      }
      else if (flag == 0){
         // fft->xo[index] = -15.0;
          fft->xo[index1] = -1.0;// * fft->xi[0][na0];
          fft->xo[index2] = -1.0;// * fft->xi[0][na0];
          fft->xo[index3] = -1.0;// * fft->xi[0][na0];
          fft->firstpass[j][0]=1;
          // fft->picked[index1] = 1.0;
      }
   /* else {
	  if (fft->picked[index0]==1.0){
		fft->xo[index]=0.0;
          }
     } */
      // fft->events(index);
    }
    }
    if (z > N3/2){
          if ((lxs+i) == 0 && j == y && (k > int(N3*7.0/8.0))  /*&& fft->firstpass[j][1]==0 */)
       //|| (j >= 49.9 && j<= 50.1)||(j >= 14.9 && j<= 15.1))
       {
     /* if(me == 0){
        if (logfile) fprintf(logfile,"Meet %d, %d, %d\n", i,j,k);
        if (screen) fprintf(screen,"Meet %d, %d, %d\n", i,j,k);
      }*/
      if (flag == 1) {
        fft->xo[index1] = 1.0;
        fft->xo[index2] = 1.0;
        fft->xo[index3] = 1.0;
     /*   if(me == 0){
          if (logfile) fprintf(logfile,"xo %lf\n", fft->xo[index]);
          if (screen) fprintf(screen,"xo %lf\n", fft->xo[index]);
        }*/
        // fft->picked[index1] = dt_kmc;
      }
      else if (flag == 0){
         // fft->xo[index] = -15.0;
          fft->xo[index1] = -1.0;// * fft->xi[0][na0];
          fft->xo[index2] = -1.0;// * fft->xi[0][na0];
          fft->xo[index3] = -1.0;// * fft->xi[0][na0];
          fft->firstpass[j][1]=1;
          // fft->picked[index1] = 1.0;
      }
   /* else {
	  if (fft->picked[index0]==1.0){
		fft->xo[index]=0.0;
          }
     } */
      // fft->events(index);
    }
    }
  }

}
/* ----------------------------------------------------------------------
   KMC events condition check
   ------------------------------------------------------------------------- */

int AppFCC3TestKMC::event_check(int i)
{


}
/* ----------------------------------------------------------------------
   Calculate KMC event time
   ------------------------------------------------------------------------- */
// 
void AppFCC3TestKMC::calculate_event_time(int index,int y)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int ND = fft->dimension;
  double * Hv = fft->material->Hv; // Pull these from the read_activation_file
  double * Pv = fft->material->Pv;  
  //int i, j, k; //added to access spatial information
  double H;
  int ida, idb;
  double kappa = 8.6173E-5; // Boltzmann constant
  beta = temperature * kappa;
  for(int i=0;i<NE;i++){

    // IMPLEMENT SOME SPATIAL PARAMETRIZATION OF ACTIVATION ENERGY HERE

    // H = Hv[y];  
    // pre_fac[i] = Pv[y];
    
    // Spatial dependence: binary regions of lo/hi activation energy 

    // if ((y < int(N2/2)+10 && y > int(N2/2)-10)) {
    //     pre_fac[i] = 6.339517E14;
    //     H = 0.2054;
    // } 
    // else {
    //     pre_fac[i] = 6.339517E14;
    //     H = 1.2054;
    // }

    // pre_fac[i] = 6.339517E14;
    // H = 0.2054;
    // rate[i] = pre_fac[i] * exp(-H/beta);


    rate[i] = Pv[y]*exp(-Hv[y]/beta);

 /* if(me == 0){
    if (logfile) fprintf(logfile,"Rate %.0e\n", rate[i]);
    if (screen) fprintf(screen,"Rate %.0e\n", rate[i]);
  }*/

  } // i

  double rannum = ranapp->uniform();
  // dt_kmc = 1 / rate[0] / 2.48 / 1.0E-14 * log(1/rannum);

  if(me == 0){
    if (logfile) fprintf(logfile,"pre_fac %.0e H %.0e\n", pre_fac[0],H);
    if (screen) fprintf(screen,"pre_fac %.0e H %.0e\n",  pre_fac[0],H);
  }
  // rate -> event time

  // compare event time, give kflag and dt_kmc
  // kflag = 1;
  //dt_kmc = 20.0;
/*  if(me == 0){
    if (logfile) fprintf(logfile,"Ran uni() %lf\n", rannum);
    if (screen) fprintf(screen,"Ran uni() %lf\n", rannum);
  }*/
  
  if(me == 0){
    // if (logfile) fprintf(logfile,"Cal_eve_time %lf\n", dt_kmc);
    // if (screen) fprintf(screen,"Cal_eve_time %lf\n", dt_kmc);
  }
}
/* ----------------------------------------------------------------------
   computes the xi in a given point and slip system
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_mean_xi(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip system, real value
    slip = 0;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 2){ // 1 slip system, imaginary value
    slip = 0;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 3){ // 2 slip system, real value
    slip = 1;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 4){ // 2 slip system, imaginary value
    slip = 1;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 5){ // 3 slip system, real value
    slip = 2;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 6){ // 3 slip system, imaginary value
    slip = 2;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  if(flag == 7){ // 4 slip system, real value
    slip = 3;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 8){ // 4 slip system, imaginary value
    slip = 3;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 9){ // 5 slip system, real value
    slip = 4;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 10){ // 5 slip system, imaginary value
    slip = 4;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 11){ // 6 slip system, real value
    slip = 5;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 12){ // 6 slip system, imaginary value
    slip = 5;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  
  return fft->xi[0][na];
}
/* ----------------------------------------------------------------------
   computes the disregistry in a given point and slip system
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_mean_dreg(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NP = fft->num_planes;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip plane, real value
    slip = 0;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 2){ // 1 slip plane, imaginary value
    slip = 0;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 3){ // 2 slip plane, real value
    slip = 1;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 4){ // 2 slip plane, imaginary value
    slip = 1;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  
  return fft->dreg[0][na];
}
/* ----------------------------------------------------------------------
   computes the stress in a given point
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_stress(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int ND = dimension;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1)  // pxx
    na = 2*(i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + 0*lN1*N2*N3*ND);
  else if(flag == 2)  // pyy
    na = 2*(i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + 1*lN1*N2*N3*ND);
  else if(flag == 3)  // pzz
    na = 2*(i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + 2*lN1*N2*N3*ND);
  else if(flag == 4)  // pxy
    na = 2*(i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + 1*lN1*N2*N3*ND);
  else if(flag == 5)  // pxz
    na = 2*(i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + 2*lN1*N2*N3*ND);
  else if(flag == 6)  // pxz
    na = 2*(i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + 2*lN1*N2*N3*ND);

  return fft->data_sigma[na];
}

/* ----------------------------------------------------------------------
   computes theta in a given point and slip system
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_theta(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }

  return fft->theta_d[na];
}
/* ----------------------------------------------------------------------
   computes nhn in a given point and slip system
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_nhn(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NP = fft->num_planes;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 2){ // 2 slip system, real value
    slip = 1;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }

  return fft->nhn_d[na]; // change for one order parameter case
}
/* ----------------------------------------------------------------------
   computes picked in a given point and slip system
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_picked(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  if(flag == 2){ // 2 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }

  return fft->picked_d[na]; //also change this for 3 order parameter case
}
/* ----------------------------------------------------------------------
   computes gradx in a given point and slip system
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_gradx(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  if(flag == 2){ // 2 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }

  return fft->gradx_d[na];
}
// /* ----------------------------------------------------------------------
//    computes gradx in a given point and slip system
//    ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_grady(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  if(flag == 2){ // 2 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }

  return fft->grady_d[na];
}
// /* ----------------------------------------------------------------------
//    computes gradx in a given point and slip system
//    ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_gradz(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  if(flag == 2){ // 2 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }

  return fft->gradz_d[na];
}
/* ----------------------------------------------------------------------
   returns local USFE
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_usfe(int n, int flag)
{
  int N1 = fft->nx;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  slip = flag - 1;
  na = (lxs+i)*N2*N3 + j*N3 + k + slip*N1*N2*N3;

  return fft->material->A_mpea[na];
}

/* ----------------------------------------------------------------------
   computes the strain in a given point
   ------------------------------------------------------------------------- */

double AppFCC3TestKMC::compute_strain(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int ND = dimension;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1)  // exx
    na = 2*(i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + 0*lN1*N2*N3*ND);
  else if(flag == 2)  // eyy
    na = 2*(i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + 1*lN1*N2*N3*ND);
  else if(flag == 3)  // ezz
    na = 2*(i*N2*N3 + j*N3 + k + 2*lN1*N2*N3 + 2*lN1*N2*N3*ND);
  else if(flag == 4)  // exy
    na = 2*(i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + 1*lN1*N2*N3*ND);
  else if(flag == 5)  // exz
    na = 2*(i*N2*N3 + j*N3 + k + 0*lN1*N2*N3 + 2*lN1*N2*N3*ND);
  else if(flag == 6)  // exz
    na = 2*(i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + 2*lN1*N2*N3*ND);

  return fft->data_eps[na];
}

/* ----------------------------------------------------------------------
   computes the xo in a given point and slip system
   ------------------------------------------------------------------------- */
double AppFCC3TestKMC::compute_xo(int n, int flag)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int NS = fft->slip_systems;
  int slip = 0;
  int na=-1;

  int i=siteijk[n][0] - lxs;
  int j=siteijk[n][1];
  int k=siteijk[n][2];

  if(flag == 1){ // 1 slip system, real value
    slip = 0;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 2){ // 2 slip system, real value
    slip = 1;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 3){ // 3 slip system, real value
    slip = 2;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  if(flag == 4){ // 4 slip system, real value
    slip = 3;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 5){ // 5 slip system, real value
    slip = 4;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 6){ // 6 slip system, real value
    slip = 5;
    na = (i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  
  return fft->xo[na];
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppFCC3TestKMC::stats(char *strtmp)
{
  char *strpnt = strtmp;
  sprintf(strpnt," %10g %10d",time,steps);
  strpnt += strlen(strpnt);

  // for (int m = 0; m < nspecies; m++) {
  //   sprintf(strpnt," %d",pcount[m]);
  //   strpnt += strlen(strpnt);
  // }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppFCC3TestKMC::stats_header(char *strtmp)
{
  char *strpnt = strtmp;
  sprintf(strpnt," %10s %10s","Time","Step");
  strpnt += strlen(strpnt);

  // for (int m = 0; m < nspecies; m++) {
  //   sprintf(strpnt," %s",sname[m]);
  //   strpnt += strlen(strpnt);
  // }
}
