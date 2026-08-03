/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_GL_X_ObsFree.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "timer.h"
#include "fft.h"
#include "dfftw_mpi.h"
#include "app.h"
#include "output.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

SolveGLXObsFree::SolveGLXObsFree(PFDD_C *pfdd_p, int narg, char **arg) :
Solve(pfdd_p, narg, arg)
{
  if (narg < 1) error->all(FLERR,"Illegal solve command");

  random = new RandomPark(ranmaster->uniform());
  it = 0;
  max_iter = 10;
  sstate = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"max_iter") == 0) {
      iarg += 1;
      max_iter = atoi(arg[iarg]);
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"tol") == 0) {
      iarg += 1;
      tol = atof(arg[iarg]);
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"aflag") == 0) {
      iarg += 1;
      aflag = atoi(arg[iarg]);
      iarg += 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

SolveGLXObsFree::~SolveGLXObsFree()
{
  delete random;
  //memory->destroy(prob);
}

/* ---------------------------------------------------------------------- */


void SolveGLXObsFree::iterate()
{
  int ND = app->dimension;
  int NS = fft->slip_systems;
  int stopsteps = app->stopsteps;
  int vol = fft->nx*fft->ny*fft->nz;
  int eff_vol = fft->eff_vol;
  timer->stamp();

  if(me == 0){
    printf("Beginning Time Evolution\n");
  }

  for(sstate = 0; sstate<stopsteps; sstate++){  // Loop on stress states
    fft->rotate_stress();
    app->resolSS();

    // GL minimization
    // if(sstate > 0) max_iter = 1;
    for(int it=0; it<max_iter; it++){
      // Initialize fft arrays
      fft->init_loop();
      // Gradient and theta calculation
      if(app->dflag == 1){
        fft->gradient_fcc();
      }
      else{
        fft->gradient();
      }
      // Calculate the core energy
      app->core_energy();
      //
      // Application specific fft forward
      fft->prep_forward();

      //Compute optimal xi from fft
      fft->internal_energy();

      // Take inverse FFT
      fft->prep_backward();
      // Update Order Parameter
      fft->update_order_parameter();

      fft->xslip_check();

      if (fft->xinorm < tol || fft->ximin > 1 || it == max_iter-1){
        if(me==0){
          printf("printing output it %d\n",it);
        }
        app->nextoutput = output->compute(app->time,1);
      }
      else if (app->time >= app->nextoutput){
        app->nextoutput = output->compute(app->time,0);
      }
      timer->stamp(TIME_OUTPUT);

      fft->prepare_next_itr();

      //compute delta, ddelta
      // fft->project_core_energy();

      timer->stamp(TIME_FFT);

      app->time++;


      if(me==0){
      printf("xinorm %f\n",fft->xinorm);
      }

      // if(me==0)
      //   printf("ximin %f\n",fft->ximin);


      // Exit if converged
      if(fft->xinorm < tol){
        if(me==0)
          printf("Loading step %d converged after %d iterations\n",sstate,it);
        break;
      }
      // if(fft->ximin > 1){
      //   if(me==0)
      //     printf("Minimum order parameter > 1 on slip plane (exiting)\n");
      //   it=max_iter;
      //   sstate=stopsteps;
      // }
      if(fft->Hflag >= 1.0 && sstate % 2 != 0){
        if(me == 0) printf("CROSSED BARRIER, EXITING\n");
        exit(1);
      }
      if((fft->X1flag >= 1.0 || fft->X2flag >= 1.0) && sstate % 2 != 0){
      for (int is=0; is<NS; is++){
       fft->tau[is] = 0.0;
        // fft->sigma[i][j] = fft->sigma[i][j] + fft->deltasig[i][j];
    }
        for (int i=0;i<ND;i++){
      for (int j=0;j<ND;j++){
        fft->sigma[i][j] = 0.0;
        fft->sigma_rot[i][j] = 0.0;
      }
    }
        // if(me == 0) printf("CROSS SLIP DETECTED, EXITING\n");
        // exit(1);
      }    
    }
      // if(fft->xinorm >= tol){
      //   if(me==0)
      //     fprintf(stderr,"Exiting step %d DID NOT CONVERGE\n",sstate);
      //   exit(1);
      // }
    if(fft->Hflag >= 0.0 && fft->X1flag >= 0.0 && fft->X2flag >= 0.0 && sstate % 2 == 0){
      fft->screw_center();
      fft->xslip_sxtal_6SFCCPlanar();
    }
    // // if(fft->xinorm >= 10.0) exit(1);
    if(fft->Hflag >= 0.0 && fft->X1flag >= 0.0 && fft->X2flag >= 0.0 && sstate % 2 != 0) app->initial_sxtal();
    // increment sigma
    for (int i=0; i<ND; i++){
      for (int j=0; j<ND; j++){
        if(sstate % 2 == 0) fft->sigma[i][j] = fft->sigma[i][j] + ((sstate/2)+1)*fft->deltasig[i][j];
        else fft->sigma[i][j] = 0.0;
        // fft->sigma[i][j] = fft->sigma[i][j] + fft->deltasig[i][j];
      }
    }
  }
  timer->stamp(TIME_SOLVE);
}
