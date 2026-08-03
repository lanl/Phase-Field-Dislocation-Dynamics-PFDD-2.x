/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef DIAG_CLASS
DiagStyle(strain,DiagStrain)
#else

#ifndef PFDD_DIAG_STRAIN_H
#define PFDD_DIAG_STRAIN_H

#include "stdio.h"
#include "diag.h"
#include <vector>

using namespace std;

namespace PFDD_NS {

  class DiagStrain : public Diag {
  public:
    DiagStrain(class PFDD_C *, int, char **);
    ~DiagStrain() {}
    void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:

  double **gstrain;          // Global stress and strain

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
