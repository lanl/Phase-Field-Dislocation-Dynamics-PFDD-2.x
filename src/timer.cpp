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

#include "mpi.h"
#include "timer.h"
#include "memory.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

Timer::Timer(PFDD_C *pfdd_p) : Pointers(pfdd_p)
{
  array = (double *) memory->smalloc(TIME_N*sizeof(double),"array");
}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
  memory->sfree(array);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
  for (int i = 0; i < TIME_N; i++) array[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void Timer::stamp()
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  previous_time = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp(int which)
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  double current_time = MPI_Wtime();
  array[which] += current_time - previous_time;
  previous_time = current_time;
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start(int which)
{
  MPI_Barrier(world);
  array[which] = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop(int which)
{
  MPI_Barrier(world);
  double current_time = MPI_Wtime();
  array[which] = current_time - array[which];
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(int which)
{
  double current_time = MPI_Wtime();
  return (current_time - array[which]);
}
