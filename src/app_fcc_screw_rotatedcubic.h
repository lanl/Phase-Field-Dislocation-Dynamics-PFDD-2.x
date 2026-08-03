/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(FCCScrewRotatedCubic,AppFCCScrewRotatedCubic)

#else

#ifndef PFDD_APP_FCCScrewRotatedCubic
#define PFDD_APP_FCCScrewRotatedCubic

#include "app.h"

namespace PFDD_NS {

  class AppFCCScrewRotatedCubic : public App {
  public:
    AppFCCScrewRotatedCubic(class PFDD_C *, int, char **);
    virtual ~AppFCCScrewRotatedCubic();
    void grow_app();
    void init_app();
    void input_app(char *, int, char **);
    void set_slip();
    double compute_mean_xi(int,int);
    double compute_mean_intE(int,int);
    double compute_data_real(int,int);
    double compute_stress(int,int);
    double compute_astress(int,int);
    double compute_tstress(int,int);
    double compute_pstress(int,int);
    double compute_rstress(int,int);
    double compute_strain(int,int);
    double compute_istrain(int,int);
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
