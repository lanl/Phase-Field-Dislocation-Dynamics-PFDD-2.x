/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(3d-hcp-basal,App3DhcpBasal)

#else

#ifndef PFDD_APP_3D_HCP_B_H
#define PFDD_APP_3D_HCP_B_H

#include "app.h"

namespace PFDD_NS {

  class App3DhcpBasal : public App {
  public:
    App3DhcpBasal(class PFDD_C *, int, char **);
    virtual ~App3DhcpBasal();
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
    double oflag;
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
