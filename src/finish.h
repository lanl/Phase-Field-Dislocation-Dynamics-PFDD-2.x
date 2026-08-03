/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifndef PFDD_FINISH_H
#define PFDD_FINISH_H

#include "pointers.h"

namespace PFDD_NS {

class Finish : protected Pointers {
 public:
  Finish(class PFDD_C *, int);
  ~Finish() {}

 private:
  void stats(int, double *, double *, double *, double *, int, int *);
};

}

#endif
