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

#ifdef APP_CLASS
AppStyle(1L1S110Bfcc,App1L1S110Bfcc)

#else

#ifndef PFDD_APP_1L1S110Bfcc_H
#define PFDD_APP_1L1S110Bfcc_H

#include "app.h"

namespace PFDD_NS {

  class App1L1S110Bfcc : public App {
  public:
    App1L1S110Bfcc(class PFDD_C *, int, char **);
    virtual ~App1L1S110Bfcc();
    void grow_app();
    void init_app();
    void input_app(char *, int, char **);
    void set_slip();
    double compute_mean_xi(int,int);
    double compute_stress(int,int);
    double compute_strain(int,int);
    void stats(char *);
    void stats_header(char *);

  protected:
    int oflag;
  };

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

E: Invalid site index
Site index must be positive and smaller or equal than the total
number of sites.

*/
