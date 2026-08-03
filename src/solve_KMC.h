/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef SOLVE_CLASS
SolveStyle(KMC,SolveKMC)

#else

#ifndef PFDD_SOLVE_KMC_H
#define PFDD_SOLVE_KMC_H

#include "solve.h"

namespace PFDD_NS {

class SolveKMC : public Solve {
 public:

  double sum;           // Sum of propensities
  int num_active;       // number of sites active
  
  SolveKMC(class PFDD_C *, int, char **);
  ~SolveKMC();
  SolveKMC *clone();

  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

  void iterate();

// This was private. Changed to protected. Inheritance (solve_linear_null)
 protected:
  class RandomPark *random;
  int nevents;
  double *prob;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

*/
