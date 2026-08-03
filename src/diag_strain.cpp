/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#include "mpi.h"
#include "diag_strain.h"
#include "error.h"
#include "app.h"
#include "comm_lattice.h"
#include "string.h"
#include <vector>
#include "fft.h"
#include "output.h"

using namespace std;
using namespace PFDD_NS;

DiagStrain::DiagStrain(PFDD_C *pfdd_p, int narg, char **arg) :
  Diag(pfdd_p,narg,arg) {

  gstrain = NULL;
}

void DiagStrain::init()
{

}

/* ----------------------------------------------------------------------
   Calculates the global strain
 ------------------------------------------------------------------------- */

void DiagStrain::compute()
{
  int ND = fft->dimension;
  if(!output->comp_strain){
    for(int i=0;i<ND;i++){
      for(int j=0;j<ND;j++){
        fft->avepst[i][j] = 0.0;
        fft->avepsts[i][j] = 0.0;
        fft->avepsd[i][j] = 0.0;
        fft->aveps[i][j] = 0.0;
        fft->ave_eps[i][j] = 0.0;
        fft->aveps_NoExt[i][j] = 0.0;
        fft->ave_eps_NoExt[i][j] = 0.0;
        if (app->nphase>1) {fft->aveVs[i][j] = 0.0;}
      }
    }
    fft->strain();
    fft->total_average_strain();
    output->comp_strain = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DiagStrain::stats(char *strtmp)
{
  sprintf(strtmp," %8g %8g %8g %8g %8g %8g %8g %8g %8g %8g %8g %8g",fft->ave_eps[0][0],fft->ave_eps[1][1],fft->ave_eps[2][2],fft->ave_eps[0][1],fft->ave_eps[0][2],fft->ave_eps[1][2],
  fft->ave_eps_NoExt[0][0],fft->ave_eps_NoExt[1][1],fft->ave_eps_NoExt[2][2],fft->ave_eps_NoExt[0][1],fft->ave_eps_NoExt[0][2],fft->ave_eps_NoExt[1][2]);
}

/* ---------------------------------------------------------------------- */

void DiagStrain::stats_header(char *strtmp)
{
  sprintf(strtmp," %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s","e11","e22","e33","e12","e13","e23","e11NE","e22NE","e33NE","e12NE","e13NE","e23NE");
}
