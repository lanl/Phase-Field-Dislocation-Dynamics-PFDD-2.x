/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef SOLVE_CLASS
SolveStyle(GL,SolveGL)

#else

#ifndef PFDD_SOLVE_GL_H
#define PFDD_SOLVE_GL_H

#include "solve.h"

namespace PFDD_NS {

class SolveGL : public Solve {
 public:
  SolveGL(class PFDD_C *, int, char **);
  ~SolveGL();
  
  void iterate();

  // This was private. Changed to protected. Inheritance (solve_linear_null)
 protected:
  class RandomPark *random;
  
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
