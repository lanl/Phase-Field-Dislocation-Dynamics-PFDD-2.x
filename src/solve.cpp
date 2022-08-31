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


  
