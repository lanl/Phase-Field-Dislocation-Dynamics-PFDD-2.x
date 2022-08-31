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

  Solve(class PFDD_C *, int, char **);
  virtual ~Solve();

  // pure virtual functions, must be defined in child class
  virtual void iterate() = 0;

};

}

#endif
