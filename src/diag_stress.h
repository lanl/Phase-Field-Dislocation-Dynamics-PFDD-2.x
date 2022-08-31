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

#ifdef DIAG_CLASS
DiagStyle(stress,DiagStress)
#else

#ifndef PFDD_DIAG_STRESS_H
#define PFDD_DIAG_STRESS_H

#include "stdio.h"
#include "diag.h"
#include <vector>

using namespace std;

namespace PFDD_NS {

  class DiagStress : public Diag {
  public:
    DiagStress(class PFDD_C *, int, char **);
    ~DiagStress() {}
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
