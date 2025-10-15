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

#ifdef SOLVE_CLASS
SolveStyle(GLCond,SolveGLCond)

#else

#ifndef PFDD_SOLVE_GLCond_H
#define PFDD_SOLVE_GLCond_H

#include "solve.h"

namespace PFDD_NS {

class SolveGLCond : public Solve {
 public:
  SolveGLCond(class PFDD_C *, int, char **);
  ~SolveGLCond();
  
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
