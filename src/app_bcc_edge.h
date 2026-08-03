/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(BCCEdge,AppBCCEdge)

#else

#ifndef PFDD_APP_BCCEdge
#define PFDD_APP_BCCEdge

#include "app.h"

namespace PFDD_NS {

  class AppBCCEdge : public App {
  public:
    AppBCCEdge(class PFDD_C *, int, char **);
    virtual ~AppBCCEdge();
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
