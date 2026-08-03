/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#include "mpi.h"
#include "diag_energy.h"
#include "error.h"
#include "app.h"
#include "comm_lattice.h"
#include "string.h"
#include <vector>
#include "fft.h"
#include "output.h"

using namespace std;
using namespace PFDD_NS;

DiagEnergy::DiagEnergy(PFDD_C *pfdd_p, int narg, char **arg) :
  Diag(pfdd_p,narg,arg) {

}

void DiagEnergy::init()
{

}

/* ----------------------------------------------------------------------
   Calculates the global strain
 ------------------------------------------------------------------------- */

void DiagEnergy::compute()
{
  int nsize = fft->nx*fft->ny*fft->nz;
  int eff_vol = fft->eff_vol;
  fft->intave = 0;
  fft->extave = 0;
  fft->latave = 0;
  if(output->comp_energy == 0){
    fft->energy_calc();
    output->comp_energy = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats(char *strtmp)
{
  sprintf(strtmp," %8g %8g %8g",fft->intave,fft->extave,fft->latave);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats_header(char *strtmp)
{
  sprintf(strtmp," %8s %8s %8s","Int","Ext","Lat");
}
