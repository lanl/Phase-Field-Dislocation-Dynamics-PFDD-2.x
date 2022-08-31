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

#ifndef PFDD_DIAG_H
#define PFDD_DIAG_H

#include "pointers.h"

namespace PFDD_NS {

class Diag : protected Pointers {
 public:
  char *style;
  int stats_flag;                   // 1 if stats drives output, 0 if not
  double next_time,delta;           // output params for stats_flag = 0
  double scale,delay;
  int logfreq,nrepeat;

  Diag(class PFDD_C *, int, char **);
  virtual ~Diag();

  // pure virtual functions, must be defined in child class
  
  virtual void init() = 0;
  virtual void compute() = 0;

  // virtual functions, may be overridden in child class

  virtual void stats(char *strtmp) {strtmp[0] = '\0';};
  virtual void stats_header(char *strtmp) {strtmp[0] = '\0';};

 protected:
  int me,nprocs;
  int iarg_child;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

*/
