/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifndef PFDD_RANDOM_MARS_H
#define PFDD_RANDOM_MARS_H

#include "pointers.h"

namespace PFDD_NS {

class RanMars : protected Pointers {
 public:
  RanMars(class PFDD_C *);
  ~RanMars();
  void init(int);
  double uniform();

 private:
  int initflag;
  int i97,j97;
  double c,cd,cm;
  double *u;
};

}

#endif

/* ERROR/WARNING messages:

E: Seed command has not been used

The seed command must be used if another command requires random
numbers.

*/
