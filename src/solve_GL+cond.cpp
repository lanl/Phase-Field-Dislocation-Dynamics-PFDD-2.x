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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_GL+cond.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "timer.h"
#include "fft.h"
#include "app.h"
#include "output.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

SolveGLCond::SolveGLCond(PFDD_C *pfdd_p, int narg, char **arg) :
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

SolveGLCond::~SolveGLCond()
{
  delete random;
  //memory->destroy(prob);
}

/* ---------------------------------------------------------------------- */


void SolveGLCond::iterate()
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
      fft->gradient();
      // Calculate the core energy
      app->core_energy();
      //
      // Application specific fft forward
      fft->prep_forward();
      fft->forwardcond();

      //Compute optimal xi from fft
      fft->internal_energy();
      
      //Compute divergence of electric current
      fft->divcurrent();

      // Take inverse FFT
      fft->prep_backward();
      fft->backwardcond();
      printf("\nconductivity_local: %f \n", double(nprocs)*fft->conductivity);
      if (me == 0) printf("\nconductivity: %f \n", fft->conductivity_global);

      // Update Order Parameters
      fft->update_order_parameter();
      fft->update_conduct_order_param();

      if (app->time >= app->nextoutput)
	app->nextoutput = output->compute(app->time,0);
      timer->stamp(TIME_OUTPUT);

      fft->prepare_next_itr();
      fft->prepare_next_conditr();

      //compute delta, ddelta
      if(app->non_Schmid)
      	fft->project_core_energy();

      timer->stamp(TIME_FFT);

      app->time++;

      // Exit if converge
      if(fft->xinorm < tol && fft->condxinorm < 1.e-13*tol)
      fft->prep_backward();
      fft->backwardcond();
      printf("\nconductivity_local: %f \n", double(nprocs)*fft->conductivity);
      if (me == 0) printf("\nconductivity: %f \n", fft->conductivity_global);
	break;
      if(fft->condxinorm > 1e-2 || fft->conductivity < 0.){
        printf("\nERROR: condxi not converging!\n");
        break;}
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
