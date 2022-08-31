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
#include "create_box.h"
#include "app.h"
#include "fft.h"
#include "region.h"
#include "error.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(PFDD_C *pfdd_p) : Pointers(pfdd_p) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
  if (fft == NULL) error->all(FLERR,"Create_box command before fft_style set");

  if (narg != 1) error->all(FLERR,"Illegal create_box command");

  if (fft->box_exist) 
    error->all(FLERR,"Cannot create box after simulation box is defined");
  if (fft->dimension == 2 && fft->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");
  if (fft->dimension == 1 && 
      (fft->yperiodic == 0 || fft->zperiodic == 0))
    error->all(FLERR,"Cannot run 1d simulation with nonperiodic Y or Z dimension");

  // region check

  int iregion = fft->find_region(arg[0]);
  if (iregion == -1) error->all(FLERR,"Create_box region ID does not exist");
  if (fft->regions[iregion]->interior == 0)
    error->all(FLERR,"Create_box region must be of type inside");

  // setup simulation box from region extent

  fft->boxxlo = fft->regions[iregion]->extent_xlo;
  fft->boxxhi = fft->regions[iregion]->extent_xhi;
  fft->boxylo = fft->regions[iregion]->extent_ylo;
  fft->boxyhi = fft->regions[iregion]->extent_yhi;
  fft->boxzlo = fft->regions[iregion]->extent_zlo;
  fft->boxzhi = fft->regions[iregion]->extent_zhi;

  fft->set_box();
  fft->box_exist = 1;

  if (fft->me == 0)
    if (screen) {
      if (screen) fprintf(screen,"Created box = (%g %g %g) to (%g %g %g)\n",
			  fft->boxxlo,fft->boxylo,fft->boxzlo,
			  fft->boxxhi,fft->boxyhi,fft->boxzhi);
      if (logfile) fprintf(logfile,"Created box = (%g %g %g) to (%g %g %g)\n",
			   fft->boxxlo,fft->boxylo,fft->boxzlo,
			   fft->boxxhi,fft->boxyhi,fft->boxzhi);
    }

  if (fft->dimension == 1) fft->procs2domain_1d();
  if (fft->dimension == 2) fft->procs2domain_2d();
  if (fft->dimension == 3) fft->procs2domain_3d();
  
  if (fft->me == 0)
    if (screen) {
      if (screen) fprintf(screen,"  %d by %d by %d processor grid\n",
			  fft->procgrid[0],fft->procgrid[1],
			  fft->procgrid[2]);
      if (logfile) fprintf(logfile,"  %d by %d by %d processor grid\n",
			   fft->procgrid[0],fft->procgrid[1],
			   fft->procgrid[2]);
    }
}
