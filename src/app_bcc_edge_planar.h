/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(BCCEdgePlanar,AppBCCEdgePlanar)

#else

#ifndef PFDD_APP_BCCEdgePlanar
#define PFDD_APP_BCCEdgePlanar

#include "app.h"

namespace PFDD_NS {

  class AppBCCEdgePlanar: public App {
  public:
    AppBCCEdgePlanar(class PFDD_C *, int, char **);
    virtual ~AppBCCEdgePlanar();
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
