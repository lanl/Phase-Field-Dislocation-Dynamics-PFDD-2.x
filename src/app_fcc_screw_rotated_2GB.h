/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(FCCScrewRotated2GB,AppFCCScrewRotated2GB)

#else

#ifndef PFDD_APP_FCCScrewRotated2GB
#define PFDD_APP_FCCScrewRotated2GB

#include "app.h"

namespace PFDD_NS {

  class AppFCCScrewRotated2GB : public App {
  public:
    AppFCCScrewRotated2GB(class PFDD_C *, int, char **);
    virtual ~AppFCCScrewRotated2GB();
    void grow_app();
    void init_app();
    void input_app(char *, int, char **);
    void set_slip();
    double compute_mean_xi(int,int);
    double compute_stress(int,int);
    double compute_astress(int,int);
    double compute_tstress(int,int);
    double compute_pstress(int,int);
    double compute_rstress(int,int);
    double compute_strain(int,int);
    double compute_istrain(int,int);
    double compute_theta(int,int);
    double compute_usfe(int, int);
    double compute_phase(int,int);
    double compute_delta(int,int);
    double compute_ddelta(int,int);
    void stats(char *);
    void stats_header(char *);

  protected:
    int oflag;
  };

}

#endif
#endif
