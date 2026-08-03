/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef SOLVE_CLASS
SolveStyle(GLXObsFree,SolveGLXObsFree)

#else

#ifndef PFDD_SOLVE_GL_X_OBSFREE_H
#define PFDD_SOLVE_GL_X_OBSFREE_H

#include "solve.h"

namespace PFDD_NS {

class SolveGLXObsFree : public Solve {
 public:
  SolveGLXObsFree(class PFDD_C *, int, char **);
  ~SolveGLXObsFree();
  
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
