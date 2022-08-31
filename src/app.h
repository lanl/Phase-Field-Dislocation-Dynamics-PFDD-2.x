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
  double stoptime;        // time at which to stop this run
  int steps;              // # of steps to run
  int stopsteps;              // # of steps to run
  int sites_exist;        // 1 if sites have been created
  double temperature;     // Temperature
  double t_inverse;       // 1/kT
  double timestep;        // timestep
  int core_flag;        // Flag for the core energy model: 1 perfect, 2 extended, 3
  int setup_flag;         //Flag for initial dislocation setup
  int dimension;        // Dimensionality of the calculation
  double nextoutput;
  int non_Schmid;                      // null for Schmid calculation (default)
  double angle_to_110, w1,w2,w3;        //MRSSP angle, non-Schmid coefficient

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
  virtual double compute_stress(int, int) {return 0;}
  virtual double compute_strain(int, int) {return 0;}
  virtual double compute_theta(int, int) {return 0;}
  virtual double compute_usfe(int, int) {return 0;}
  virtual double compute_delta(int, int) {return 0;}
  virtual double compute_ddelta(int, int) {return 0;}

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
