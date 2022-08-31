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

#ifdef COMMAND_CLASS
CommandStyle(set,Set)

#else

#ifndef PFDD_SET_H
#define PFDD_SET_H

#include "pointers.h"

namespace PFDD_NS {

class Set : protected Pointers {
 public:
  Set(class PFDD_C *);
  void command(int, char **);

 private:
  int siteindex;
  int count;
  int ivalue,ivaluelo,ivaluehi;
  double dvalue,dvaluelo,dvaluehi;
  int loopflag,regionflag,iregion;
  double fraction;

  struct Condition {                     // list of if-test conditions
    int lhs,type,index,stride;
    int op;
    int irhs;
    double drhs;
  };
  Condition *cond;
  int ncondition;

  int latticeflag;
  class App *app;

    
  void set_single(int, int);
  void set_range(int, int);
  void set_displace(int, int);
  int condition(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Set command before sites exist

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Setting a quantity application does not support

The application defines what variables it supports.  You cannot set a
variable with the set command on a variable that isn't supported.

E: Set command region ID does not exist

Self-explanatory.

E: Set if test on quantity application does not support

The application defines what variables it supports.  You cannot do an
if test with the set command on a variable that isn't supported.

*/
