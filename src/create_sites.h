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

#ifdef COMMAND_CLASS
CommandStyle(create_sites,CreateSites)

#else

#ifndef PFDD_CREATE_SITES_H
#define PFDD_CREATE_SITES_H

#include "pointers.h"

namespace PFDD_NS {

class CreateSites : protected Pointers {
  friend class ReadSites;

 public:
  CreateSites(class PFDD_C *);
  void command(int, char **);
  void ghosts_from_connectivity(class App *, int);

 private:
  int style,nregion,valueflag,valueindex,ivalue;
  double dvalue;
  int *basisflag,*basis_ivalue;
  double *basis_dvalue;
  int maxneigh;

  int latticeflag;
  //class App *app;
  
  tagint **idneigh;            // global indices of neighbors of each site
                               // same as AppLattice neighbor, but tagint,
                               // tmp storage until convert to local indices

  int **siteijk;               // global indices of each site
                               // 0,1,2 = i,j,k lattice indices
                               // 3 = which basis atom in unit cell

  int ***cmap;                 // connectivity map for regular lattices
                               // cmap[nbasis][maxneigh][4]
                               // 0,1,2 = i,j,k lattice unit cell offsets
                               // 3 = which basis atom in unit cell


  // geometry info for building structured lattice with neighbors

  int nx,ny,nz;
  int xlo,xhi,ylo,yhi,zlo,zhi;
  double xorig,yorig,zorig;
  int latstyle,nbasis;
  double xlattice,ylattice,zlattice;

  struct Site {
    int id,proc,index;
    double x,y,z;
  };

  void structured_lattice();
  void structured_connectivity();

  void offsets(double **);
  void offsets_2d(int, double **, double, double, int, int **);
  void offsets_3d(int, double **, double, double, int, int **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Create_sites command before app_style set

Self-explanatory.

E: Create_sites command before simulation box is defined

Self-explanatory.

E: Cannot create sites after sites already exist

Self-explanatory.

E: Cannot create sites with undefined lattice

Must use lattice commands first to define a lattice.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Create_sites region ID does not exist

Self-explanatory.

E: Creating a quantity application does not support

The application defines what variables it supports.  You cannot set a
variable with the create_sites command for a variable that isn't
supported.

E: Must use value option before basis option in create_sites command

Self-explanatory.

E: Cannot use create_sites basis with random lattice

Self-explanatory.

E: Periodic box is not a multiple of lattice spacing

UNDOCUMENTED

E: Did not create correct number of sites

One or more created sites were not assigned to a processor
correctly.

E: Bad neighbor site ID

UNDOCUMENTED

E: Random lattice has no connectivity

The cutoff distance is likely too short.

E: Ghost site was not found

Internal SPPARKS error.  Should not occur.

E: Ghost connection was not found

Internal SPPARKS error.  Should not occur.

E: Incorrect lattice neighbor count

Internal SPPARKS error.

*/
