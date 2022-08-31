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
