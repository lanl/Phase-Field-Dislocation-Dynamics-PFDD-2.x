/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifndef PFDD_REGION_H
#define PFDD_REGION_H

#include "pointers.h"

namespace PFDD_NS {

class Region : protected Pointers {
 public:
  char *id,*style;
  int interior;                     // 1 for interior, 0 for exterior
  double xscale,yscale,zscale;      // scale factors for lattice units
  double extent_xlo,extent_xhi;     // bounding box on region
  double extent_ylo,extent_yhi;
  double extent_zlo,extent_zhi;
  
  Region(class PFDD_C *, int, char **);
  virtual ~Region();
  virtual int match(double, double, double) = 0;

 protected:
  void options(int, char **);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Use of region with undefined lattice

The lattice command must be used before defining a geometric region.

*/
