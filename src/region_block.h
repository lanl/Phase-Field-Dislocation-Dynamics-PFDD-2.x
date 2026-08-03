/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef REGION_CLASS
RegionStyle(block,RegBlock)

#else

#ifndef PFDD_REGION_BLOCK_H
#define PFDD_REGION_BLOCK_H

#include "region.h"

namespace PFDD_NS {

class RegBlock : public Region {
 public:
  RegBlock(class PFDD_C *, int, char **);
  int match(double, double, double);

 private:
  double xlo,xhi,ylo,yhi,zlo,zhi;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use region INF or EDGE when box does not exist

Can only define a region with these parameters after a simulation
box has been defined.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

*/
