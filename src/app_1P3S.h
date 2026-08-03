/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(1P3S,App1P3S)

#else

#ifndef PFDD_APP_1P3S
#define PFDD_APP_1P3S

#include "app.h"

namespace PFDD_NS {

  class App1P3S : public App {
  public:
    App1P3S(class PFDD_C *, int, char **);
    virtual ~App1P3S();
    void grow_app();
    void init_app();
    void input_app(char *, int, char **);
    void set_slip();
    double compute_mean_xi(int,int);
    double compute_stress(int,int);
    double compute_strain(int,int);
    double compute_theta(int,int);
    double compute_usfe(int, int);
    void stats(char *);
    void stats_header(char *);

  protected:
    int oflag;
  };

}

#endif
#endif
