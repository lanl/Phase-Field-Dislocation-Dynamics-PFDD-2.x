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

#include "mpi.h"
#include "string.h"
#include "PFDD.h"
#include "input.h"
#include "memory.h"
#include "error.h"
#include "Types.h"
#include "timer.h"
#include "random_mars.h"
#include "output.h"
#include "fft.h"
#include "app.h"
#include "solve.h"

using namespace PFDD_NS;

/* ----------------------------------------------------------------------
   start up LAMMPS
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
------------------------------------------------------------------------- */

PFDD_C::PFDD_C(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory(this);
  error = new Error(this);
  
  world = communicator;
  screen = stdout;
  logfile = NULL;

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int partscreenflag = 0;
  int partlogflag = 0;
  int cudaflag = -1;
  int helpflag = 0;
  suffix = NULL;
  suffix_enable = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-in") == 0 ||
               strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-screen") == 0 ||
               strcmp(arg[iarg],"-sc") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-log") == 0 ||
               strcmp(arg[iarg],"-l") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-var") == 0 ||
               strcmp(arg[iarg],"-v") == 0) {
      if (iarg+3 > narg)
        error->all(FLERR,"Invalid command-line argument");
      iarg += 3;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
    } else if (strcmp(arg[iarg],"-echo") == 0 ||
               strcmp(arg[iarg],"-e") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-pscreen") == 0 ||
               strcmp(arg[iarg],"-ps") == 0) {
      if (iarg+2 > narg)
       error->all(FLERR,"Invalid command-line argument");
      partscreenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-plog") == 0 ||
               strcmp(arg[iarg],"-pl") == 0) {
      if (iarg+2 > narg)
       error->all(FLERR,"Invalid command-line argument");
      partlogflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-cuda") == 0 ||
               strcmp(arg[iarg],"-c") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Invalid command-line argument");
      if (strcmp(arg[iarg+1],"on") == 0) cudaflag = 1;
      else if (strcmp(arg[iarg+1],"off") == 0) cudaflag = 0;
      else error->all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-suffix") == 0 ||
               strcmp(arg[iarg],"-sf") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Invalid command-line argument");
      delete [] suffix;
      int n = strlen(arg[iarg+1]) + 1;
      suffix = new char[n];
      strcpy(suffix,arg[iarg+1]);
      suffix_enable = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-help") == 0 ||
               strcmp(arg[iarg],"-h") == 0) {
      if (iarg+1 > narg)
        error->all(FLERR,"Invalid command-line argument");
      helpflag = 1;
      iarg += 1;
    } else error->all(FLERR,"Invalid command-line argument");
  }

    // set universe screen and logfile
  MPI_Comm_rank(world,&me);
  
  if (me == 0) {
    if (screenflag == 0)
      screen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      screen = NULL;
    else {
      screen = fopen(arg[screenflag],"w");
      if (screen == NULL)
        error->one(FLERR,"Cannot open universe screen file");
    }
    if (logflag == 0) {
      logfile = fopen("log.pfdd","w");
      if (logfile == NULL)
        error->one(FLERR,"Cannot open log.pfdd");
    } else if (strcmp(arg[logflag],"none") == 0)
      logfile = NULL;
    else {
      logfile = fopen(arg[logflag],"w");
      if (logfile == NULL)
        error->one(FLERR,"Cannot open universe log file");
    }
  }

  if (me > 0) {
    if (screenflag == 0) screen = stdout;
    else screen = NULL;
    logfile = NULL;
  }

  // make universe and single world the same, since no partition switch
  // world inherits settings from universe
  // set world screen, logfile, communicator, infile
  // open input script if from file

  infile = NULL;
  
  if (me == 0) {
    if (inflag == 0) infile = stdin;
    else infile = fopen(arg[inflag],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[inflag]);
      error->one(FLERR,str);
    }
  }
    
  // allocate input class now that MPI is fully setup

  input = new Input(this,narg,arg);

  // allocate top-level classes

  create();
  post_create();

  // if helpflag set, print help and quit

  // if (helpflag) {
  //   if (universe->me == 0) print_styles();
  //   error->done();
  // }
}

/* ----------------------------------------------------------------------
   shutdown LAMMPS
   delete top-level classes
   close screen and log files in world and universe
   output files were already closed in destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

PFDD_C::~PFDD_C()
{
  destroy();

  if (logfile) fclose(logfile);
  if (screen && screen != stdout) fclose(screen);
  
  delete input;
  delete error;
  delete memory;
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in constructor
   some classes have package variants
------------------------------------------------------------------------- */

void PFDD_C::create()
{
  // Comm class must be created before Atom class
  // so that nthreads is defined when create_avec invokes grow()

  app = NULL;
  solve = NULL;
  fft = NULL;
  ranmaster = new RanMars(this);
  output = new Output(this);
  timer = new Timer(this);


//   if (cuda) comm = new CommCuda(this);
//   else comm = new Comm(this);

//   if (cuda) neighbor = new NeighborCuda(this);
//   else neighbor = new Neighbor(this);

//   if (cuda) domain = new DomainCuda(this);
// #ifdef LMP_USER_OMP
//   else domain = new DomainOMP(this);
// #else
//   else domain = new Domain(this);
// #endif

//   atom = new Atom(this);
//   atom->create_avec("atomic",0,NULL,suffix);

//   group = new Group(this);
//   force = new Force(this);    // must be after group, to create temperature

//   if (cuda) modify = new ModifyCuda(this);
//   else modify = new Modify(this);

//   output = new Output(this);  // must be after group, so "all" exists
//                               // must be after modify so can create Computes
//   update = new Update(this);  // must be after output, force, neighbor
//   timer = new Timer(this);
}

/* ----------------------------------------------------------------------
   invoke package-specific setup commands
   called from PFDD_C constructor and after clear() command
   only invoke if suffix is set and enabled
------------------------------------------------------------------------- */

void PFDD_C::post_create()
{
  if (suffix && suffix_enable) {
    if (strcmp(suffix,"gpu") == 0) input->one("package gpu force/neigh 0 0 1");
    if (strcmp(suffix,"omp") == 0) input->one("package omp *");
  }
}


/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void PFDD_C::destroy()
{

  delete app;
  delete solve;
  delete fft;
  delete ranmaster;
  delete output;
  delete timer;

}

/* ----------------------------------------------------------------------
   for each style, print name of all child classes build into executable
------------------------------------------------------------------------- */

void PFDD_C::print_styles()
{
//   printf("\nList of style options included in this executable:\n\n");

//   printf("Atom styles:");
// #define ATOM_CLASS
// #define AtomStyle(key,Class) printf(" %s",#key);
// #include "style_atom.h"
// #undef ATOM_CLASS
//   printf("\n\n");

//   printf("Integrate styles:");
// #define INTEGRATE_CLASS
// #define IntegrateStyle(key,Class) printf(" %s",#key);
// #include "style_integrate.h"
// #undef INTEGRATE_CLASS
//   printf("\n\n");

//   printf("Minimize styles:");
// #define MINIMIZE_CLASS
// #define MinimizeStyle(key,Class) printf(" %s",#key);
// #include "style_minimize.h"
// #undef MINIMIZE_CLASS
//   printf("\n\n");

//   printf("Pair styles:");
// #define PAIR_CLASS
// #define PairStyle(key,Class) printf(" %s",#key);
// #include "style_pair.h"
// #undef PAIR_CLASS
//   printf("\n\n");

//   printf("Bond styles:");
// #define BOND_CLASS
// #define BondStyle(key,Class) printf(" %s",#key);
// #include "style_bond.h"
// #undef BOND_CLASS
//   printf("\n\n");

//   printf("Angle styles:");
// #define ANGLE_CLASS
// #define AngleStyle(key,Class) printf(" %s",#key);
// #include "style_angle.h"
// #undef ANGLE_CLASS
//   printf("\n\n");

//   printf("Dihedral styles:");
// #define DIHEDRAL_CLASS
// #define DihedralStyle(key,Class) printf(" %s",#key);
// #include "style_dihedral.h"
// #undef DIHEDRAL_CLASS
//   printf("\n\n");

//   printf("Improper styles:");
// #define IMPROPER_CLASS
// #define ImproperStyle(key,Class) printf(" %s",#key);
// #include "style_improper.h"
// #undef IMPROPER_CLASS
//   printf("\n\n");

//   printf("KSpace styles:");
// #define KSPACE_CLASS
// #define KSpaceStyle(key,Class) printf(" %s",#key);
// #include "style_kspace.h"
// #undef KSPACE_CLASS
//   printf("\n\n");

//   printf("Fix styles (upper case are only for internal use):");
// #define FIX_CLASS
// #define FixStyle(key,Class) printf(" %s",#key);
// #include "style_fix.h"
// #undef FIX_CLASS
//   printf("\n\n");

//   printf("Compute styles:");
// #define COMPUTE_CLASS
// #define ComputeStyle(key,Class) printf(" %s",#key);
// #include "style_compute.h"
// #undef COMPUTE_CLASS
//   printf("\n\n");

//   printf("Region styles:");
// #define REGION_CLASS
// #define RegionStyle(key,Class) printf(" %s",#key);
// #include "style_region.h"
// #undef REGION_CLASS
//   printf("\n\n");

//   printf("Dump styles:");
// #define DUMP_CLASS
// #define DumpStyle(key,Class) printf(" %s",#key);
// #include "style_dump.h"
// #undef DUMP_CLASS
//   printf("\n\n");

//   printf("Command styles (add-on input script commands):");
// #define COMMAND_CLASS
// #define CommandStyle(key,Class) printf(" %s",#key);
// #include "style_command.h"
// #undef COMMAND_CLASS
//   printf("\n");
}
