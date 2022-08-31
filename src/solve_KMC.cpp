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
#include "solve_KMC.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

SolveKMC::SolveKMC(PFDD_C *pfdd_p, int narg, char **arg) : 
  Solve(pfdd_p, narg, arg)
{
  if (narg != 1) error->all(FLERR,"Illegal solve command");

  random = new RandomPark(ranmaster->uniform());
  prob = NULL;
}

/* ---------------------------------------------------------------------- */

SolveKMC::~SolveKMC()
{
  delete random;
  memory->destroy(prob);
}

/* ---------------------------------------------------------------------- */

SolveKMC *SolveKMC::clone()
{
  int narg = 1;
  char *arg[1];
  arg[0] = style;

  SolveKMC *ptr = new SolveKMC(pfdd_p,narg,arg);

  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveKMC::init(int n, double *propensity)
{
  delete [] prob;
  nevents = n;
  memory->create(prob,n,"solve/linear:prob");

  sum = 0.0;
  num_active = 0;

  for (int i = 0; i < n; i++) {
    if (propensity[i] > 0.0) num_active++;
    prob[i] = propensity[i];
    sum += propensity[i];
  }
}

/* ---------------------------------------------------------------------- */

void SolveKMC::update(int n, int *indices, double *propensity)
{
  int m;
  for (int i = 0; i < n; i++) {
    m = indices[i];
    if (prob[m] > 0.0) num_active--;
    if (propensity[m] > 0.0) num_active++;
    sum -= prob[m];
    prob[m] = propensity[m];
    sum += propensity[m];
  }
}
/* ---------------------------------------------------------------------- */

void SolveKMC::update(int n, double *propensity)
{
  if (prob[n] > 0.0) num_active--;
  if (propensity[n] > 0.0) num_active++;
  sum -= prob[n];
  prob[n] = propensity[n];
  sum += propensity[n];
}
/* ---------------------------------------------------------------------- */


void SolveKMC::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}
/* ---------------------------------------------------------------------- */

int SolveKMC::event(double *pdt)
{
  int m;

  if (num_active == 0) {
    sum = 0.0;
    return -1;
  }

  double fraction = sum * random->uniform();
  double partial = 0.0;

  for (m = 0; m < nevents; m++) {
    partial += prob[m];
    if (partial > fraction) break;
  }

  *pdt = -1.0/sum * log(random->uniform());
  return m;
}

void SolveKMC::iterate()
{

}
