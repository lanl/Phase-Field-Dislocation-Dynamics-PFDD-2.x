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

#include "Types.h"
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "lattice.h"
#include "fft.h"
#include "memory.h"
#include "error.h"

using namespace PFDD_NS;

// same as in create_sites.cpp and diag_cluster.cpp

enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
     FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D,BCC_OCTA_TETRA,BCC_OCTA};

/* ---------------------------------------------------------------------- */

Lattice::Lattice(PFDD_C *pfdd_p, int narg, char **arg) : Pointers(pfdd_p)
{
  // parse style arg

  if (narg < 1) error->all(FLERR,"Illegal lattice command");

  if (strcmp(arg[0],"none") == 0) style = NONE;
  else if (strcmp(arg[0],"line/2n") == 0) style = LINE_2N;
  else if (strcmp(arg[0],"sq/4n") == 0) style = SQ_4N;
  else if (strcmp(arg[0],"sq/8n") == 0) style = SQ_8N;
  else if (strcmp(arg[0],"tri") == 0) style = TRI;
  else if (strcmp(arg[0],"sc/6n") == 0) style = SC_6N;
  else if (strcmp(arg[0],"sc/26n") == 0) style = SC_26N;
  else if (strcmp(arg[0],"fcc") == 0) style = FCC;
  else if (strcmp(arg[0],"bcc") == 0) style = BCC;
  else if (strcmp(arg[0],"diamond") == 0) style = DIAMOND;
  else if (strcmp(arg[0],"fcc/octa/tetra") == 0) style = FCC_OCTA_TETRA;
  else if (strcmp(arg[0],"bcc/octa/tetra") == 0) style = BCC_OCTA_TETRA;
  else if (strcmp(arg[0],"bcc/octa") == 0) style = BCC_OCTA;
  else if (strcmp(arg[0],"random/1d") == 0) style = RANDOM_1D;
  else if (strcmp(arg[0],"random/2d") == 0) style = RANDOM_2D;
  else if (strcmp(arg[0],"random/3d") == 0) style = RANDOM_3D;
  else error->all(FLERR,"Illegal lattice command");

  if (style == NONE) {
    if (narg > 1) error->all(FLERR,"Illegal lattice command");
    return;
  }

  if (style == LINE_2N || style == SQ_4N || style == SQ_8N ||
      style == TRI || style == SC_6N || style == SC_26N ||
      style == FCC || style == BCC || style == DIAMOND || 
      style == FCC_OCTA_TETRA || style == BCC_OCTA_TETRA ||
      style == BCC_OCTA) {
    if (narg != 2) error->all(FLERR,"Illegal lattice command");
    latconst = atof(arg[1]);
  }

  if (style == RANDOM_1D || style == RANDOM_2D || style == RANDOM_3D) {
    if (narg != 3) error->all(FLERR,"Illegal lattice command");
    latconst = 1.0;
    nrandom = ATOTAGINT(arg[1]);
    cutoff = atof(arg[2]);
  }

  // check dimensionality

  if ((style == LINE_2N || style == RANDOM_1D) && 
      fft->dimension != 1)
    error->all(FLERR,"Lattice style does not match dimension");
  if ((style == SQ_4N || style == SQ_8N || style == TRI || 
       style == RANDOM_2D) && 
      fft->dimension != 2)
    error->all(FLERR,"Lattice style does not match dimension");
  if ((style == SC_6N || style == SC_26N || style == FCC || 
       style == BCC || style == DIAMOND || style == FCC_OCTA_TETRA || style == BCC_OCTA_TETRA || style == BCC_OCTA ||
       style == RANDOM_3D) && 
      fft->dimension != 3)
    error->all(FLERR,"Lattice style does not match dimension");

  // set basis atoms for each style

  nbasis = 0;
  basis = NULL;

  if (style == LINE_2N || style == SQ_4N || style == SQ_8N ||
      style == SC_6N || style == SC_26N) {
    add_basis(0.0,0.0,0.0);
  } else if (style == TRI) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
  } else if (style == BCC) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.5);
  } else if (style == FCC) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.0,0.5,0.5);
    add_basis(0.5,0.0,0.5);
    add_basis(0.5,0.5,0.0);
  } else if (style == DIAMOND) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.0,0.5,0.5);
    add_basis(0.5,0.0,0.5);
    add_basis(0.5,0.5,0.0);
    add_basis(0.25,0.25,0.25);
    add_basis(0.25,0.75,0.75);
    add_basis(0.75,0.25,0.75);
    add_basis(0.75,0.75,0.25);
  } else if (style == FCC_OCTA_TETRA) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.0,0.5,0.5);
    add_basis(0.5,0.0,0.5);
    add_basis(0.5,0.5,0.0);
    add_basis(0.5,0.0,0.0);
    add_basis(0.0,0.5,0.0);
    add_basis(0.0,0.0,0.5);
    add_basis(0.5,0.5,0.5);
    add_basis(0.25,0.25,0.25);
    add_basis(0.75,0.25,0.25);
    add_basis(0.25,0.75,0.25);
    add_basis(0.75,0.75,0.25);
    add_basis(0.25,0.25,0.75);
    add_basis(0.75,0.25,0.75);
    add_basis(0.25,0.75,0.75);
    add_basis(0.75,0.75,0.75);
  }
  else if (style == BCC_OCTA_TETRA) {
    add_basis(0.0,0.0,0.0);   //lattice
    add_basis(0.5,0.5,0.5);
    add_basis(0.0,0.5,0.5);   //octahedral
    add_basis(0.5,0.0,0.5);
    add_basis(0.5,0.5,0.0);
    add_basis(0.5,0.0,0.0);
    add_basis(0.0,0.5,0.0);
    add_basis(0.0,0.0,0.5);
    add_basis(0.5,0.25,0.0);  //tetrahedral
    add_basis(0.75,0.5,0.0);
    add_basis(0.5,0.75,0.0);
    add_basis(0.25,0.5,0.0);
    add_basis(0.5,0.0,0.25);
    add_basis(0.75,0.0,0.5);
    add_basis(0.5,0.0,0.75);
    add_basis(0.25,0.0,0.5);
    add_basis(0.0,0.5,0.25);
    add_basis(0.0,0.75,0.5);
    add_basis(0.0,0.5,0.75);
    add_basis(0.0,0.25,0.5);
  }
  else if (style == BCC_OCTA) {
    add_basis(0.0,0.0,0.0);   //lattice
    add_basis(0.5,0.5,0.5);
    add_basis(0.0,0.5,0.5);   //octahedral
    add_basis(0.5,0.0,0.5);
    add_basis(0.5,0.5,0.0);
    add_basis(0.5,0.0,0.0);
    add_basis(0.0,0.5,0.0);
    add_basis(0.0,0.0,0.5);
  }

  // set defaults for optional args

  origin[0] = origin[1] = origin[2] = 0.0;

  orientx[0] = 1;  orientx[1] = 0;  orientx[2] = 0;
  orienty[0] = 0;  orienty[1] = 1;  orienty[2] = 0;
  orientz[0] = 0;  orientz[1] = 0;  orientz[2] = 1;

  a1[0] = 1.0;  a1[1] = 0.0;  a1[2] = 0.0;
  a2[0] = 0.0;  a2[1] = 1.0;  a2[2] = 0.0;
  a3[0] = 0.0;  a3[1] = 0.0;  a3[2] = 1.0;

  if (style == TRI) a2[1] = sqrt(3.0);

  // lattice spacings

  xlattice = a1[0]*latconst;
  ylattice = a2[1]*latconst;
  zlattice = a3[2]*latconst;
}

/* ---------------------------------------------------------------------- */

Lattice::~Lattice()
{
  memory->destroy(basis);
}

/* ----------------------------------------------------------------------
   add a basis atom to list
   x,y,z = fractional coords within unit cell
------------------------------------------------------------------------- */

void Lattice::add_basis(double x, double y, double z)
{
  memory->grow(basis,nbasis+1,3,"lattice:basis");
  basis[nbasis][0] = x;
  basis[nbasis][1] = y;
  basis[nbasis][2] = z;
  nbasis++;
}

