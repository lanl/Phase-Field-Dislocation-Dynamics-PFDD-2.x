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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "diag.h"
#include "app.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

Diag::Diag(PFDD_C *pfdd_p, int narg, char **arg) : Pointers(pfdd_p)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  // stats_flag = 1, provide stats, compute interval controlled by Output
  // stats_flag = 0, do not provide stats to Output

  stats_flag = 1;
  delta = 0.0;
  logfreq = 0;
  delay = 0.0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stats") == 0) {
      iarg++;
      if (iarg < narg) {
	if (strcmp(arg[iarg],"yes") == 0) stats_flag = 1;
	else if (strcmp(arg[iarg],"no") == 0) stats_flag = 0;
	else error->all(FLERR,"Illegal diag_style command");
      } else error->all(FLERR,"Illegal diag_style command");
    } else if (strcmp(arg[iarg],"delta") == 0) {
      iarg++;
      if (iarg < narg) {
	delta = atof(arg[iarg]);
	if (delta <= 0.0) error->all(FLERR,"Illegal diag_style command");
      } else error->all(FLERR,"Illegal diag_style command");
    } else if (strcmp(arg[iarg],"logfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal diag_style command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (scale <= 0) error->all(FLERR,"Illegal diag_style command");
      if (nrepeat < 0) error->all(FLERR,"Illegal diag_style command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"loglinfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal diag_style command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (scale <= 0) error->all(FLERR,"Illegal diag_style command");
      if (nrepeat < 0) error->all(FLERR,"Illegal diag_style command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 2;
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      iarg++;
      if (iarg < narg) {
	delay = atof(arg[iarg]);
      } else error->all(FLERR,"Illegal diag_style command");
    } else break;
    iarg++;
  }

  iarg_child = iarg;

  if (stats_flag && logfreq) error->all(FLERR,"Illegal diag_style command");
  if (stats_flag && delta > 0.0) error->all(FLERR,"Illegal diag_style command");
  if (stats_flag && delay > 0.0) error->all(FLERR,"Illegal diag_style command");
}

/* ---------------------------------------------------------------------- */

Diag::~Diag()
{
  delete [] style;
}
