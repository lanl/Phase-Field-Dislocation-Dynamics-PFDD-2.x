/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#include "string.h"
#include "solve.h"
#include "fft.h"
#include <math.h>
#include "mpi.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

Solve::Solve(PFDD_C *pfdd_p, int narg, char **arg) : Pointers(pfdd_p)
{
  int n = strlen(arg[0]) + 1;
  style = NULL;
  style = new char[n];
  strcpy(style,arg[0]);

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

}

/* ---------------------------------------------------------------------- */

Solve::~Solve()
{
  delete [] style;
}

/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */


  
