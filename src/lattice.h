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

#ifndef PFDD_LATTICE_H
#define PFDD_LATTICE_H

#include "pointers.h"

namespace PFDD_NS {

class Lattice : protected Pointers {
 public:
  int style;                           // enum list of NONE,SC,FCC,etc
  double xlattice,ylattice,zlattice;   // lattice scale factors in 3 dims
  double a1[3],a2[3],a3[3];            // edge vectors of unit cell
  int nbasis;                          // # of basis atoms in unit cell
  double **basis;                      // fractional coords of each basis atom
                                       // within unit cell (0 <= coord < 1)
  tagint nrandom;                      // # of sites for random lattices
  double cutoff;                       // neighbor cutoff for random lattices

  Lattice(class PFDD_C *, int, char **);
  ~Lattice();

private:
  double latconst;                     // lattice constant
  double origin[3];                    // lattice origin
  int orientx[3];                      // lattice orientation vecs
  int orienty[3];                      // orientx = what lattice dir lies
  int orientz[3];                      //           along x dim in box

  void add_basis(double, double, double);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Lattice style does not match dimension

Self-explanatory.

E: Cannot use coloring without domain nx,ny,nz defined

UNDOCUMENTED

E: Color stencil is incommensurate with lattice size

Since coloring induces a pattern of colors, this pattern
must fit an integer number of times into a periodic lattice.

*/
