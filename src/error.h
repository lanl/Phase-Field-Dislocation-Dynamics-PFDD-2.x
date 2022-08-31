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

#ifndef PFDD_ERROR_H
#define PFDD_ERROR_H

#include "pointers.h"

namespace PFDD_NS {

class Error : protected Pointers {
 public:
  Error(class PFDD_C *);

  void all(const char *, int, const char *);
  void one(const char *, int, const char *);
  void warning(const char *, int, const char *, int = 1);
  void message(const char *, int, const char *, int = 1);
  void done();
};

}

#endif
