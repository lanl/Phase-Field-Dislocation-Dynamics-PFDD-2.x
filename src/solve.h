/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifndef PFDD_SOLVE_H
#define PFDD_SOLVE_H

#include "pointers.h"

namespace PFDD_NS {

class Solve : protected Pointers {
 public:
  char *style;
  int me, nprocs;

  int max_iter;         // Max # of iteration in the minimization procedure
  int it;               // iteration in the minimization (time)
  int sstate;           // iteration in the stress
  double tol;           // Tolerance for inner loop convergence
  int aflag;

  Solve(class PFDD_C *, int, char **);
  virtual ~Solve();

  // pure virtual functions, must be defined in child class
  virtual void iterate() = 0;

};

}

#endif
