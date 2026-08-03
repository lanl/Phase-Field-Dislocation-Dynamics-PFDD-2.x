/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
   ------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_3SCTBTest.h"
#include "fft.h"
#include "solve.h"
#include "random_park.h"
#include "error.h"

using namespace PFDD_NS;;

/* ---------------------------------------------------------------------- */

App3SCTBTest::App3SCTBTest(PFDD_C *pfdd_p, int narg, char **arg) :
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

  if(setup_flag == -1)
    setup_flag = 3; //default to a screw dipole
  if(core_flag == -1)
    core_flag = 1; //default to sin^2
}

/* ---------------------------------------------------------------------- */

App3SCTBTest::~App3SCTBTest()
{
  delete [] sites;
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
   ------------------------------------------------------------------------- */

void App3SCTBTest::grow_app()
{
  darray = fft->xi;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
   ------------------------------------------------------------------------- */

void App3SCTBTest::init_app()
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
void App3SCTBTest::input_app(char *command, int narg, char **arg)
{
  // error->all(FLERR,"Unrecognized command");
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


void App3SCTBTest::set_slip()
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
// MUST BE IN THIS ORDER FOR LOOKUP TABLE OR FOURIER: <110> <101> <011>

  fft->xb[0][0] = -1.0/s2;
  fft->xb[0][1] = 1.0/s2;
  fft->xb[0][2] = 0.0;

  fft->xb[1][0] = 1.0/s2;
  fft->xb[1][1] = 0.0;
  fft->xb[1][2] = -1.0/s2;

  fft->xb[2][0] = 0.0;
  fft->xb[2][1] = -1.0/s2;
  fft->xb[2][2] = 1.0/s2;

}
/* ----------------------------------------------------------------------
   computes the xi in a given point and slip system
   ------------------------------------------------------------------------- */

double App3SCTBTest::compute_mean_xi(int n, int flag)
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

double App3SCTBTest::compute_mean_dreg(int n, int flag)
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
   computes the internal energy in a given point and slip system
------------------------------------------------------------------------- */

double App3SCTBTest::compute_mean_intE(int n, int flag)
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
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  else if(flag == 2){ // 2 slip system, real value
    slip = 1;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  else if(flag == 3){ // 3 slip system, real value
    slip = 2;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  if(flag == 4){ // 4 slip system, real value
    slip = 3;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  else if(flag == 5){ // 5 slip system, real value
    slip = 4;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  else if(flag == 6){ // 6 slip system, real value
    slip = 5;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }

  return fft->E_int[na];
}
/* ----------------------------------------------------------------------
   computes the external energy in a given point and slip system
------------------------------------------------------------------------- */

double App3SCTBTest::compute_mean_latE(int n, int flag)
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
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  else if(flag == 2){ // 2 slip system, real value
    slip = 1;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  else if(flag == 3){ // 3 slip system, real value
    slip = 2;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  if(flag == 4){ // 4 slip system, real value
    slip = 3;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  else if(flag == 5){ // 5 slip system, real value
    slip = 4;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  else if(flag == 6){ // 6 slip system, real value
    slip = 5;
    na = i*N2*N3 + j*N3 + k + slip*lN1*N2*N3;;
  }
  return fft->E_lat[na];
}
/* ----------------------------------------------------------------------
   computes the lattice energy in a given point and slip system
------------------------------------------------------------------------- */

double App3SCTBTest::compute_mean_extE(int n, int flag)
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
  }
  else if(flag == 2){ // 2 slip system, real value
    slip = 1;
  }
  else if(flag == 3){ // 3 slip system, real value
    slip = 2;
  }
  if(flag == 4){ // 4 slip system, real value
    slip = 3;
  }
  else if(flag == 5){ // 5 slip system, real value
    slip = 4;
  }
  else if(flag == 6){ // 6 slip system, real value
    slip = 5;
  }

  return fft->E_ext[na];
}
/* ----------------------------------------------------------------------
   computes the stress in a given point
   ------------------------------------------------------------------------- */

double App3SCTBTest::compute_stress(int n, int flag)
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
  else if(flag == 6)  // pyz
    na = 2*(i*N2*N3 + j*N3 + k + 1*lN1*N2*N3 + 2*lN1*N2*N3*ND);

  return fft->data_sigma[na];
}

/* ----------------------------------------------------------------------
   computes theta in a given point and slip system
   ------------------------------------------------------------------------- */

double App3SCTBTest::compute_theta(int n, int flag)
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
   returns local USFE
   ------------------------------------------------------------------------- */

double App3SCTBTest::compute_usfe(int n, int flag)
{
  int N1 = fft->nx;
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

  // slip = flag - 1;
  // na = (lxs+i)*N2*N3 + j*N3 + k + slip*N1*N2*N3;

  // return fft->material->A_mpea[na];

  na = (flag*N1*N2*N3 + (lxs+i)*N2*N3 + j*N3 + k + slip*7*N1*N2*N3);
  // if (me == 0) printf("fcc7_mpea: %f\n",)
  return fft->material->fcc7_mpea[na];
}

/* ----------------------------------------------------------------------
   computes the strain in a given point
   ------------------------------------------------------------------------- */

double App3SCTBTest::compute_strain(int n, int flag)
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
double App3SCTBTest::compute_xo(int n, int flag)
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
   computes nhn in a given point and slip system
   ------------------------------------------------------------------------- */

double App3SCTBTest::compute_nhn(int n, int flag)
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

double App3SCTBTest::compute_picked(int n, int flag)
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
   print stats
------------------------------------------------------------------------- */

void App3SCTBTest::stats(char *strtmp)
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

void App3SCTBTest::stats_header(char *strtmp)
{
  char *strpnt = strtmp;
  sprintf(strpnt," %10s %10s","Time","Step");
  strpnt += strlen(strpnt);

  // for (int m = 0; m < nspecies; m++) {
  //   sprintf(strpnt," %s",sname[m]);
  //   strpnt += strlen(strpnt);
  // }
}
