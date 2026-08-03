/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
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
