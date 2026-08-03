/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef DIAG_CLASS
DiagStyle(energy,DiagEnergy)
#else

#ifndef PFDD_DIAG_ENERGY_H
#define PFDD_DIAG_ENERGY_H

#include "stdio.h"
#include "diag.h"
#include <vector>

using namespace std;

namespace PFDD_NS {

  class DiagEnergy : public Diag {
  public:
    DiagEnergy(class PFDD_C *, int, char **);
    ~DiagEnergy() {}
    void init();
  void compute();
  void stats(char *);
  void stats_header(char *);


};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid diag style

Check the input script syntax and compare to the documentation
for the command.

E: Invalid spin value for magnetization

The spin of one site must be 1 or 2.

*/
