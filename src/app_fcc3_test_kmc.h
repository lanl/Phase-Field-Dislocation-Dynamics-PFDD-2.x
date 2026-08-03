/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(FCC3TestKMC,AppFCC3TestKMC)

#else

#ifndef PFDD_APP_FCC3_TEST_KMC
#define PFDD_APP_FCC3_TEST_KMC

#include "app.h"

namespace PFDD_NS {

  class AppFCC3TestKMC : public App {
  public:
    AppFCC3TestKMC(class PFDD_C *, int, char **);
    virtual ~AppFCC3TestKMC();
    void grow_app();
    void init_app();
    void input_app(char *, int, char **);
    void set_slip();
    void events();
    void kmc_stress_update(int,int);
    void barrier(int,int,int);
    int event_check(int);
    void calculate_event_time(int,int); //Changed
    double compute_mean_xi(int,int);
    double compute_mean_dreg(int,int);
    double compute_stress(int,int);
    double compute_strain(int,int);
    double compute_theta(int,int);
    double compute_usfe(int,int);
    double compute_xo(int,int);
    double compute_nhn(int,int);
    double compute_picked(int,int);
    double compute_gradx(int,int);
    double compute_grady(int,int);
    double compute_gradz(int,int);
    void stats(char *);
    void stats_header(char *);
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
