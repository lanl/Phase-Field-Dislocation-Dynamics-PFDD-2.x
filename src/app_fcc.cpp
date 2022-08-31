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
#include "string.h"
#include "stdlib.h"
#include "app_fcc.h"
#include "fft.h"
#include "solve.h"
#include "random_park.h"
#include "error.h"

using namespace PFDD_NS;;

/* ---------------------------------------------------------------------- */

AppFCC::AppFCC(PFDD_C *pfdd_p, int narg, char **arg) :
  App(pfdd_p,narg,arg)
{
  ndouble = 0;
  timestep = 1.0;
  //neighshell = 1;
  slip_systems = 12;  // The # of slip systems has to be allocated in the constructure
                      //  fot the FFT class to be able to allocate arrays
  CD = 0.5;           // Dislocation mobility

  if(setup_flag == -1)
    setup_flag = 4; //default to mode1
  if(core_flag == -1)
    core_flag = 1; //default to sin^2
}

/* ---------------------------------------------------------------------- */

AppFCC::~AppFCC()
{
  delete [] sites;
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
   ------------------------------------------------------------------------- */

void AppFCC::grow_app()
{
  darray = fft->xi;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
   ------------------------------------------------------------------------- */

void AppFCC::init_app()
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
void AppFCC::input_app(char *command, int narg, char **arg)
{
  error->all(FLERR,"Unrecognized command");
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


void AppFCC::set_slip()
{
  int i, j, k;
  double s3;
  double xnf=0,ynf=0,znf=0;
  double xbf=0,ybf=0,zbf=0;

  s3 = sqrt(3);
  for (i=0;i<3 ; i++) {
    fft->xn[i][0]= 1.0/s3;
    fft->xn[i][1]= 1.0/s3;
    fft->xn[i][2]= 1.0/s3;
  }

  for (i=3;i<6 ; i++) {
    fft->xn[i][0]= -1.0/s3;
    fft->xn[i][1]= 1.0/s3;
    fft->xn[i][2]= 1.0/s3;
  }

  for (i=6;i<9 ; i++) {
    fft->xn[i][0]= 1.0/s3;
    fft->xn[i][1]= -1.0/s3;
    fft->xn[i][2]= 1.0/s3;
  }
  for (i=9;i<12 ; i++) {
    fft->xn[i][0]= 1.0/s3;
    fft->xn[i][1]= 1.0/s3;
    fft->xn[i][2]= -1.0/s3;
  }
  fft->xb[0][0] = -1.0/2.0;
  fft->xb[0][1] = 1.0/2.0;
  fft->xb[0][2] = 0.0;

  fft->xb[1][0] = 0.0;
  fft->xb[1][1] = -1.0/2.0;
  fft->xb[1][2] = 1.0/2.0;

  fft->xb[2][0] = 1.0/2.0;
  fft->xb[2][1] = 0.0;
  fft->xb[2][2] = -1.0/2.0;

  fft->xb[3][0] = -1.0/2.0;
  fft->xb[3][1] = -1.0/2.0;
  fft->xb[3][2] = 0.0;

  fft->xb[4][0] = 1.0/2.0;
  fft->xb[4][1] = 0.0;
  fft->xb[4][2] = 1.0/2.0;

  fft->xb[5][0] = 0.0;
  fft->xb[5][1] = -1.0/2.0;
  fft->xb[5][2] = 1.0/2.0;

  fft->xb[6][0] = -1.0/2.0;
  fft->xb[6][1] = -1.0/2.0;
  fft->xb[6][2] = 0.0;

  fft->xb[7][0] = 1.0/2.0;
  fft->xb[7][1] = 0.0;
  fft->xb[7][2] = -1.0/2.0;

  fft->xb[8][0] = 0.0;
  fft->xb[8][1] = -1.0/2.0;
  fft->xb[8][2] = -1.0/2.0;

  fft->xb[9][0] = -1.0/2.0;
  fft->xb[9][1] = 1.0/2.0;
  fft->xb[9][2] = 0.0;

  fft->xb[10][0] = 1.0/2.0;
  fft->xb[10][1] = 0.0;
  fft->xb[10][2] = 1.0/2.0;

  fft->xb[11][0] = 0.0;
  fft->xb[11][1] = -1.0/2.0;
  fft->xb[11][2] = -1.0/2.0;

  // Rotation
  for(int i=0; i<slip_systems; i++){
    xnf = (fft->xn[i][0]*fft->four[0][0] + fft->xn[i][1]*fft->four[0][1] + fft->xn[i][2]*fft->four[0][2])
      /fft->normfour[0];
    ynf = (fft->xn[i][0]*fft->four[1][0] + fft->xn[i][1]*fft->four[1][1] + fft->xn[i][2]*fft->four[1][2])
      /fft->normfour[1];
    znf = (fft->xn[i][0]*fft->four[2][0] + fft->xn[i][1]*fft->four[2][1] + fft->xn[i][2]*fft->four[2][2])
      /fft->normfour[2];

    fft->xn[i][0] = xnf;
    fft->xn[i][1] = ynf;
    fft->xn[i][2] = znf;

    xbf = (fft->xb[i][0]*fft->four[0][0] + fft->xb[i][1]*fft->four[0][1] + fft->xb[i][2]*fft->four[0][2])
      /fft->normfour[0];
    ybf = (fft->xb[i][0]*fft->four[1][0] + fft->xb[i][1]*fft->four[1][1] + fft->xb[i][2]*fft->four[1][2])
      /fft->normfour[1];
    zbf = (fft->xb[i][0]*fft->four[2][0] + fft->xb[i][1]*fft->four[2][1] + fft->xb[i][2]*fft->four[2][2])
      /fft->normfour[2];

    fft->xb[i][0] = xbf;
    fft->xb[i][1] = ybf;
    fft->xb[i][2] = zbf;
  }

}

/* ----------------------------------------------------------------------
   computes the xi in a given point and slip system
   ------------------------------------------------------------------------- */

double AppFCC::compute_mean_xi(int n, int flag)
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
  else if(flag == 7){ // 4 slip system, real value
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
  else if(flag == 13){ // 7 slip system, real value
    slip = 6;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 14){ // 7 slip system, imaginary value
    slip = 6;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 15){ // 8 slip system, real value
    slip = 7;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 16){ // 8 slip system, imaginary value
    slip = 7;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 17){ // 9 slip system, real value
    slip = 8;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 18){ // 9 slip system, imaginary value
    slip = 8;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 19){ // 10 slip system, real value
    slip = 9;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 20){ // 10 slip system, imaginary value
    slip = 9;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 21){ // 11 slip system, real value
    slip = 10;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 22){ // 11 slip system, imaginary value
    slip = 10;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }
  else if(flag == 23){ // 12 slip system, real value
    slip = 11;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3);
  }
  else if(flag == 24){ // 12 slip system, imaginary value
    slip = 11;
    na = 2*(i*N2*N3 + j*N3 + k + slip*lN1*N2*N3) + 1;
  }

  return fft->xi[0][na];
}

/* ----------------------------------------------------------------------
   computes the stress in a given point
   ------------------------------------------------------------------------- */

double AppFCC::compute_stress(int n, int flag)
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
   computes the strain in a given point
   ------------------------------------------------------------------------- */

double AppFCC::compute_strain(int n, int flag)
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

void AppFCC::stats(char *strtmp)
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

void AppFCC::stats_header(char *strtmp)
{
  char *strpnt = strtmp;
  sprintf(strpnt," %10s %10s","Time","Step");
  strpnt += strlen(strpnt);

  // for (int m = 0; m < nspecies; m++) {
  //   sprintf(strpnt," %s",sname[m]);
  //   strpnt += strlen(strpnt);
  // }
}
