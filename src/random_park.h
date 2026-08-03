/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifndef PFDD_RANDOM_PARK_H
#define PFDD_RANDOM_PARK_H

#include "Types.h"

namespace PFDD_NS {

class RandomPark {
 public:
  int seed;

  RandomPark(int);
  RandomPark(double);
  ~RandomPark() {}
  void reset(double, int, int);
  void tagreset(double, tagint, int);
  double uniform();
  int irandom(int);
  tagint tagrandom(tagint);
  bigint bigrandom(bigint);
};

}

#endif
