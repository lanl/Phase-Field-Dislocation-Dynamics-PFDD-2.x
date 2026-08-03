/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
   ------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_1L1S_111B_kmc.h"
#include "fft.h"
#include "solve.h"
#include "random_park.h"
#include "error.h"
#include "memory.h"

using namespace PFDD_NS;;

/* ---------------------------------------------------------------------- */

App1L1S111BKMC::App1L1S111BKMC(PFDD_C *pfdd_p, int narg, char **arg) :
  App(pfdd_p,narg,arg)
{
  ndouble = 0;
  timestep = 1.0;
  //neighshell = 1;
  slip_systems = 1;  // The # of slip systems has to be allocated in the constructure
                      //  fot the FFT class to be able to allocate arrays
  num_planes = 1;  // The # of plane normals has to be allocated in the constructure
                      //  fot the FFT class to be able to allocate arrays

  CD = 0.5;           // Dislocation mobility
  oflag = 1;

  NE = 1;             // total # of events for KMC
  // memory->create(event_cond,NE,"event_cond");
  memory->create(event_time,NE,"event_time");
  memory->create(pre_fac,NE,"pre_fac");
  memory->create(E_act,NE,"E_act");
  memory->create(rate,NE,"rate");
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
    setup_flag = 10; //default to straight ortho
  if(core_flag == -1)
    core_flag = 1; //default to sin^2
}

/* ---------------------------------------------------------------------- */

App1L1S111BKMC::~App1L1S111BKMC()
{
  delete [] sites;
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
   ------------------------------------------------------------------------- */

void App1L1S111BKMC::grow_app()
{
  darray = fft->xi;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
   ------------------------------------------------------------------------- */

void App1L1S111BKMC::init_app()
{
  delete [] sites;
  sites = new int[1 + maxneigh];

  grow_app();
  set_slip();
  initial_sxtal();

  stopsteps = static_cast<int>(floor(stoptime/timestep));
}
/* ----------------------------------------------------------------------
   input specific for this app
   ------------------------------------------------------------------------- */
void App1L1S111BKMC::input_app(char *command, int narg, char **arg)
{

  if (strcmp(command,"events") == 0) {
    if (narg < 9) error->all(FLERR,"Illegal app_style command");
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
    if(me == 0){
      if (logfile) fprintf(logfile,"event No %d, pre_fac = %lf, E_act = %lf\n", ievent, pre_fac[ievent], E_act[ievent]);
      if (screen) fprintf(screen,"event No %d, pre_fac = %lf, E_act = %lf\n", ievent, pre_fac[ievent], E_act[ievent]);
    }
  }

}

/* ----------------------------------------------------------------------
   compute energy of site
   ------------------------------------------------------------------------- */

// double App1L1S111BKMC::site_energy(int i)
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

// void App1L1S111BKMC::site_event_rejection(int i, RandomPark *random)
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

// double App1L1S111BKMC::site_propensity(int i)
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

// void App1L1S111BKMC::site_event(int i, RandomPark *random)
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

// int App1L1S111BKMC::site_spin (int i) {
//   if ((i<0)||(i>nlocal))
//     error->all (FLERR, "Invalid site index");
//   return spin[i];
// }


void App1L1S111BKMC::set_slip()
{
  int i, j, k;
  double s3, s2;
  double xnf=0,ynf=0,znf=0;
  double xbf=0,ybf=0,zbf=0;
  s3 = sqrt(3.0);
  s2 = sqrt(2.0);

  fft->xn[0][0]= 1.0/s2;
  fft->xn[0][1]= 0;
  fft->xn[0][2]= -1.0/s2;

  fft->xb[0][0]= 1.0/s3;
  fft->xb[0][1]= 1.0/s3;
  fft->xb[0][2]= 1.0/s3;

  // Rotation
  // for(int i=0; i<slip_systems; i++){
  //   xnf = (fft->xn[i][0]*fft->four[0][0] + fft->xn[i][1]*fft->four[0][1] + fft->xn[i][2]*fft->four[0][2])
  //     /fft->normfour[0];
  //   ynf = (fft->xn[i][0]*fft->four[1][0] + fft->xn[i][1]*fft->four[1][1] + fft->xn[i][2]*fft->four[1][2])
  //     /fft->normfour[1];
  //   znf = (fft->xn[i][0]*fft->four[2][0] + fft->xn[i][1]*fft->four[2][1] + fft->xn[i][2]*fft->four[2][2])
  //     /fft->normfour[2];

  //   fft->xn[i][0] = xnf;
  //   fft->xn[i][1] = ynf;
  //   fft->xn[i][2] = znf;

  //   xbf = (fft->xb[i][0]*fft->four[0][0] + fft->xb[i][1]*fft->four[0][1] + fft->xb[i][2]*fft->four[0][2])
  //     /fft->normfour[0];
  //   ybf = (fft->xb[i][0]*fft->four[1][0] + fft->xb[i][1]*fft->four[1][1] + fft->xb[i][2]*fft->four[1][2])
  //     /fft->normfour[1];
  //   zbf = (fft->xb[i][0]*fft->four[2][0] + fft->xb[i][1]*fft->four[2][1] + fft->xb[i][2]*fft->four[2][2])
  //     /fft->normfour[2];

  //   fft->xb[i][0] = xbf;
  //   fft->xb[i][1] = ybf;
  //   fft->xb[i][2] = zbf;
  // }
  return;
}

/* ----------------------------------------------------------------------
   It creates a file with the initial
   values of the order parameters
   ------------------------------------------------------------------------- */

// void App1L1S111BKMC::initial_sxtal()
// {
//     fft->initial_sxtal_1LNonOrthoFCC();
//     // fft->initial_sxtal_StraightNonOrtho();
// }
/* ----------------------------------------------------------------------
   Check if KMC events happen and their time
   ------------------------------------------------------------------------- */

void App1L1S111BKMC::events()
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
  int na0, index0;
  for(i=0;i<lN1;i++)
  for(j=0;j<N2;j++)
  for(k=0;k<N3;k++){
    index = i*N2*N3 + j*N3 + k + is*lN1*N2*N3;
    index0 = 44*N2*N3 + j*N3 + k + is*lN1*N2*N3;
    na0 = 2 * index;
    if (fft->picked[index]>0.0 
       && (((lxs+i) >= 48.9 && (lxs+i)<= 49.1)||((lxs+i) >= 14.9 && (lxs+i)<= 15.1) )
       //|| (j >= 49.9 && j<= 50.1)||(j >= 14.9 && j<= 15.1))
       ) {
     /* if(me == 0){
        if (logfile) fprintf(logfile,"Meet %d, %d, %d\n", i,j,k);
        if (screen) fprintf(screen,"Meet %d, %d, %d\n", i,j,k);
      }*/
      if (fft->xo[index] == 0.0) {
        fft->xo[index] = 1.0;
     /*   if(me == 0){
          if (logfile) fprintf(logfile,"xo %lf\n", fft->xo[index]);
          if (screen) fprintf(screen,"xo %lf\n", fft->xo[index]);
        }*/
        calculate_event_time(index);
        // fft->picked[index] = dt_kmc;
      }
      else if (fft->xo[index] == 1.0){
        if (fft->picked[index]<=1.0) {
         // fft->xo[index] = -15.0;
          fft->xo[index] = -1.0;// * fft->xi[0][na0];
          fft->picked[index] = 1.0;
        }
        else{
          fft->picked[index] -= 1.0;
        }
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
/* ----------------------------------------------------------------------
   KMC events condition check
   ------------------------------------------------------------------------- */

int App1L1S111BKMC::event_check(int i)
{


}
/* ----------------------------------------------------------------------
   Calculate KMC event time
   ------------------------------------------------------------------------- */

void App1L1S111BKMC::calculate_event_time(int index)
{
  int lN1 = fft->local_x;
  int lxs = fft->local_x_start;
  int N2 = fft->local_y;
  int N3 = fft->local_z;
  int ND = fft->dimension;
  double H;
  int ida, idb;
  double kappa = 8.6173E-5; // Boltzmann constant
  beta = temperature * kappa;
  for(int i=0;i<NE;i++){
    /*
    H = 0.0;
    for (ida=0; ida<ND; ida++){
      for (idb=0; idb<ND; idb++){
        // add ijk loop
        index = na + 2*(ida*lN1*N2*N3 + idb*lN1*N2*N3*ND);
        H -= vol_ac[i][ida][idb] * fft->data_sigma[index];
      }
    }
    H += E_act[i];
    rate[i] = pre_fac[i] * exp(-beta*H); */

    //pre_fac[i] = 1.0E13;
   // H = 0.1;
    pre_fac[i] = 6.339517E14;
    H = 0.2054;
    //H = 1.659233;
  /*if(me == 0){
    if (logfile) fprintf(logfile,"pre_fac %.0e H %.0e\n", pre_fac[0],H);
    if (screen) fprintf(screen,"pre_fac %.0e H %.0e\n",  pre_fac[0],H);
  }*/
    rate[i] = pre_fac[i] * exp(-H/beta);

 /* if(me == 0){
    if (logfile) fprintf(logfile,"Rate %.0e\n", rate[i]);
    if (screen) fprintf(screen,"Rate %.0e\n", rate[i]);
  }*/

  } // i
  // rate -> event time

  // compare event time, give kflag and dt_kmc
  // kflag = 1;
  //dt_kmc = 20.0;
  double rannum = ranapp->uniform();
/*  if(me == 0){
    if (logfile) fprintf(logfile,"Ran uni() %lf\n", rannum);
    if (screen) fprintf(screen,"Ran uni() %lf\n", rannum);
  }*/
  // dt_kmc = 1 / rate[0] / 1.78 / 1.0E-14 * log(1/rannum);
  if(me == 0){
    // if (logfile) fprintf(logfile,"Cal_eve_time %lf\n", dt_kmc);
    // if (screen) fprintf(screen,"Cal_eve_time %lf\n", dt_kmc);
  }
}
/* ----------------------------------------------------------------------
   computes the xi in a given point and slip system
   ------------------------------------------------------------------------- */

double App1L1S111BKMC::compute_mean_xi(int n, int flag)
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

  return fft->xi[0][na];
}

/* ----------------------------------------------------------------------
   computes the stress in a given point
   ------------------------------------------------------------------------- */

double App1L1S111BKMC::compute_stress(int n, int flag)
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

double App1L1S111BKMC::compute_theta(int n, int flag)
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

  return fft->theta[na];
}
/* ----------------------------------------------------------------------
   computes nhn in a given point and slip system
   ------------------------------------------------------------------------- */

double App1L1S111BKMC::compute_nhn(int n, int flag)
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

  return fft->nhn[na];
}
/* ----------------------------------------------------------------------
   computes picked in a given point and slip system
   ------------------------------------------------------------------------- */

double App1L1S111BKMC::compute_picked(int n, int flag)
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

  return fft->picked[na];
}
/* ----------------------------------------------------------------------
   computes the strain in a given point
   ------------------------------------------------------------------------- */

double App1L1S111BKMC::compute_strain(int n, int flag)
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
   print stats
------------------------------------------------------------------------- */

void App1L1S111BKMC::stats(char *strtmp)
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

void App1L1S111BKMC::stats_header(char *strtmp)
{
  char *strpnt = strtmp;
  sprintf(strpnt," %10s %10s","Time","Step");
  strpnt += strlen(strpnt);

  // for (int m = 0; m < nspecies; m++) {
  //   sprintf(strpnt," %s",sname[m]);
  //   strpnt += strlen(strpnt);
  // }
}
