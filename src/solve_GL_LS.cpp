/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_GL_LS.h"
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

SolveGLLS::SolveGLLS(PFDD_C *pfdd_p, int narg, char **arg) :
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
  }
}

/* ---------------------------------------------------------------------- */

SolveGLLS::~SolveGLLS()
{
  delete random;
  //memory->destroy(prob);
}

/* ---------------------------------------------------------------------- */


void SolveGLLS::iterate()
{
  int ND = app->dimension;
  int stopsteps = app->stopsteps;

  timer->stamp();

  if(me == 0){
    printf("Beginning Time Evolution\n");
  }

  for(sstate = 0; sstate<stopsteps; sstate++){  // Loop on stress states
    fft->rotate_stress();
    app->resolSS();


    // GL minimization

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
      app->gb_energy();
      //
      // Application specific fft forward inverse of data_fftw
      // fft->prep_forward();

      //Compute optimal xi from fft
      fft->internal_energy();

      // Take inverse FFT of data_real
      // fft->prep_backward();

      // Update Order Parameter
      fft->update_order_parameter();

      //if (fft->ximin > 1 || it == max_iter-1){
      if (fft->xinorm < tol || fft->ximin > 1 || it == max_iter-1){
        if(me==0)
          printf("printing output it %d\n",it);
        app->nextoutput = output->compute(app->time,1);
      }
      else if (app->time >= app->nextoutput){
        app->nextoutput = output->compute(app->time,0);
      }
      timer->stamp(TIME_OUTPUT);

      fft->prepare_next_itr();

      //compute delta, ddelta
      fft->project_core_energy();

      timer->stamp(TIME_FFT);

      app->time++;


      if(me==0)
      printf("xinorm %f\n",fft->xinorm);

      if(me==0)
        printf("ximin %f\n",fft->ximin);


      // Exit if converged
      if(fft->xinorm < tol){
        if(me==0)
          printf("Loading step %d converged after %d iterations\n",sstate,it);
        break;
      }
      if(fft->ximin > 1){
        if(me==0)
          printf("Minimum order parameter > 1 on slip plane (exiting)\n");
        it=max_iter;
        sstate=stopsteps;
      }

    }
    // increment sigma
    for (int i=0; i<ND; i++){
      for (int j=0; j<ND; j++){
        fft->sigma[i][j] = fft->sigma[i][j] + fft->deltasig[i][j];
      }
    }
  }
  timer->stamp(TIME_SOLVE);
}
