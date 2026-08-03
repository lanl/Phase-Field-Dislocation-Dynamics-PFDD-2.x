/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifndef PFDD_APP_H
#define PFDD_APP_H

#include "stdio.h"
#include "pointers.h"
#include <stack>

namespace PFDD_NS {

class App : protected Pointers {
  friend class CommLattice;
 public:
  enum APP_CLASSES{GENERAL,LATTICE,OFF_LATTICE};
  int me,nprocs;

  int appclass;           // one of the enum values
  char *style;            // style name of app
  double time;            // current simulation time due to executed events
  double kmctime;
  double stoptime;        // time at which to stop this run
  int steps;              // # of steps to run
  int stopsteps;              // # of steps to run
  int sites_exist;        // 1 if sites have been created
  double temperature;     // Temperature
  double t_inverse;       // 1/kT
  double timestep;        // timestep
  int core_flag;        // Flag for the core energy model: 1 perfect, 2 extended, 3
  int setup_flag;         //Flag for initial dislocation setup
  int phase_flag;       //Flag for different phases/materials
  int dflag;            //Flag for needing disregistry evolution (FCC)
  int tflag;            //Flag for twin (Coherent Sigma3 {111} CTB)
  int xstop;            // Flag to stop stress after cross slip (SFW calc on cross slip plane)
  int dimension;        // Dimensionality of the calculation
  double nextoutput;
  int non_Schmid;                      // null for Schmid calculation (default)
  double angle_to_110, w1,w2,w3;        //MRSSP angle, non-Schmid coefficient

  int nphase;
  int gbe_flag;         //Flag for the grain boundary energy: 0 none 1 proportional to core energy
  // owned + ghost sites

  tagint nglobal;              // global # of sites
  int nlocal;                  // # of sites I own
  int nghost;                  // # of ghost sites I store
  int neighshell;              // # of nearest neighbors to consider

  int ninteger,ndouble;        // # of ints and doubles per site
  tagint *id;                  // global ID of site
  double **xyz;                // coords of site
  double **sigma;              // applied stress
  double **deltasig;           // Increment of sigma
  int **siteijk;               // global indices of each site
                               // 0,1,2 = i,j,k lattice indices
                               // 3 = which basis atom in unit cell

  int **iarray;                // one or more ints per site
  double **darray;             // one or more doubles per site
  double **xi;
  int *sites;
  double obsden;      // obstacle density
  double *xo;         // obstacles
  double CD;        // Dislocation mobility coefficient
  int slip_systems;   // # of slip systems
  int num_planes;   // # of plane normals
  double sa,sb,sc;      //scaling factors

  int nmax;                    // max # of sites per-site arrays can store
  int neighlayer;                    // cutoff layer to assign neighbors
  int maxneigh;                // max neighbors of any site in entire system
  int *numneigh;               // # of neighbors of each site
  int **neighbor;              // local indices of neighbors of each site

  // Probably I'll need a Material Class

  int NE;                       // # of events for KMC
  double beta;                  // 1/(T*kappa)
  //int *event_cond;              // if the event condition is met
  double *event_time;           // time of the event occurs
  //double cur_time;              // current time according to # of step * timestep
  int kflag;                    // flag of the event that will be performed
  double *pre_fac;
  double *E_act;
  double ***vol_ac;
  double *rate;
  double *stress_dependence; // stress dependence turned on as means of facilitation if == 1
  double dt_glide_avg,dt_glide_min,dt_xslip_avg,dt_xslip_min;

  class CommLattice *comm;
  class RandomPark *ranapp;    // RN generator for KMC and rejection KMC
  class RandomPark *ranstrict; // RN generator for per-site strict rKMC

  int *owner;                  // proc who owns the site
  int *index;                  // index of site on owning proc

  //void input(char *, int, char **);
  //void init();
  //void setup();
  //void iterate();

  void grow(int);
  void add_site(tagint, double, double, double);
  void add_ghost(tagint, double, double, double, int, int);
  void add_neighbors(int, int, char **);
  void add_values(int, char **);
  void set_initial_sxtal(char *str);
  void initial_sxtal();
  void print_connectivity();
  void bounds(char *, int, int, int &, int &);

  void input(char *, int, char **);
  void init();
  void setup();
  void iterate();

  // pure virtual functions, must be defined in child class

  virtual void grow_app() = 0;
  virtual void set_slip() = 0;

  // virtual functions, may be overridden by child class

  virtual void events() {};
  virtual void events_glide() {};
  virtual void events_xslip() {};
  virtual void events_stop() {};
  virtual void events_skip() {};
  virtual void barrier(int, int, int) {};

  virtual void input_app(char *, int, char **) {};
  virtual void init_app() {};
  virtual void setup_app() {};
  virtual void setup_end_app() {};
  virtual void iterate_app(double) {};
  virtual void *extract_app(char *) {return NULL;}

  virtual void user_update(double) {};
  virtual void core_energy();
  virtual void resolSS();
  virtual double compute_mean_xi(int, int) {return 0;}
  virtual double compute_mean_dreg(int, int) {return 0;}
  virtual double compute_mean_dregp(int, int) {return 0;}
  virtual double compute_mean_intE(int, int) {return 0;}
  virtual double compute_mean_extE(int, int) {return 0;}
  virtual double compute_mean_latE(int, int) {return 0;}
  virtual double compute_stress_NoExt(int, int) {return 0;}
  virtual double compute_data_real(int, int) {return 0;}
  virtual double compute_ff1(int, int) {return 0;}

  virtual double compute_stress(int, int) {return 0;}
  virtual double compute_astress(int, int) {return 0;}
  virtual double compute_tstress(int, int) {return 0;}
  virtual double compute_pstress(int, int) {return 0;}
  virtual double compute_ustress(int, int) {return 0;}
  virtual double compute_rstress(int, int) {return 0;}
  virtual double compute_strain(int, int) {return 0;}
  virtual double compute_theta(int, int) {return 0;}
  virtual double compute_usfe(int, int) {return 0;}
  virtual double compute_nhn(int, int) {return 0;}
  virtual double compute_picked(int, int) {return 0;}
  virtual double compute_dtkmc(int, int) {return 0;}
  virtual double compute_gradx(int, int) {return 0;}
  virtual double compute_grady(int, int) {return 0;}
  virtual double compute_gradz(int, int) {return 0;}
  virtual double compute_xo(int, int) {return 0;}
  virtual double compute_phase(int, int) {return 0;}

  virtual double compute_delta(int, int) {return 0;}
  virtual double compute_ddelta(int, int) {return 0;}
  virtual void gb_energy();

  virtual double compute_istrain(int, int) {return 0;}

  App(class PFDD_C *, int, char **);
  virtual ~App();
  void run(int, char **);
  void reset_time(double);
  void *extract(char *);
  tagint max_site_ID();

  // virtual functions, may be overridden in child class

  virtual void stats(char *strtmp) {strtmp[0] = '\0';};
  virtual void stats_header(char *strtmp) {strtmp[0] = '\0';};
  //virtual void *extract_app(char *) {return NULL;}
  virtual void set_temperature(int, char **);
  virtual void set_update_only(int, char **);
  virtual void set_sigma(int, char **);
  virtual void set_dimension(int, char **);
  virtual void set_scale(int, char **);
  virtual void set_nonSchmid(int, char **); //non-Schmid effect

 protected:
  int first_run;
  int update_only, allow_update;

  void create_arrays();
  void recreate_arrays();
  int contiguous_sites();
};

}

#endif

/* ERROR/WARNING messages:

E: Cannot run application until simulation box is defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Run upto value is before current time

Self-explanatory.

*/
