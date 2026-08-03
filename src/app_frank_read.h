/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(FrankRead,AppFrankRead)

#else

#ifndef PFDD_APP_FRANK_READ
#define PFDD_APP_FRANK_READ

#include "app.h"

namespace PFDD_NS {

  class AppFrankRead : public App {
  public:
    AppFrankRead(class PFDD_C *, int, char **);
    virtual ~AppFrankRead();
    void grow_app();
    void init_app();
    void input_app(char *, int, char **);
    void set_slip();
    void initial_sxtal();
    double compute_mean_xi(int,int);
    double compute_stress(int,int);
    double compute_strain(int,int);
    double compute_theta(int,int);
    void stats(char *);
    void stats_header(char *);

  protected:
    int oflag;
  };

}

#endif
#endif
