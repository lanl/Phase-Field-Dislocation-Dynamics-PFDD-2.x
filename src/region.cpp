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

#include "stdlib.h"
#include "string.h"
#include "region.h"
#include "fft.h"
#include "lattice.h"
#include "error.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

Region::Region(PFDD_C *pfdd_p, int narg, char **arg) : Pointers(pfdd_p)
{
  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);
}

/* ---------------------------------------------------------------------- */

Region::~Region()
{
  delete [] id;
  delete [] style;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of region input line
------------------------------------------------------------------------- */

void Region::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal region command");

  // option defaults

  interior = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal region command");
      if (strcmp(arg[iarg+1],"in") == 0) interior = 1;
      else if (strcmp(arg[iarg+1],"out") == 0) interior = 0;
      else error->all(FLERR,"Illegal region command");
      iarg += 2;
    } else error->all(FLERR,"Illegal region command");
  }

  // setup scaling

  if (fft->lattice == NULL)
    error->all(FLERR,"Use of region with undefined lattice");

  xscale = fft->lattice->xlattice;
  yscale = fft->lattice->ylattice;
  zscale = fft->lattice->zlattice;
}
