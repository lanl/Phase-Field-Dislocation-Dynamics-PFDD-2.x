/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#include "mpi.h"
#include "diag_stress.h"
#include "error.h"
#include "app.h"
#include "comm_lattice.h"
#include "string.h"
#include <vector>
#include "fft.h"
#include "output.h"

using namespace std;
using namespace PFDD_NS;

DiagStress::DiagStress(PFDD_C *pfdd_p, int narg, char **arg) :
  Diag(pfdd_p,narg,arg) {

}

void DiagStress::init()
{

}

/* ----------------------------------------------------------------------
   Calculates the global strain
 ------------------------------------------------------------------------- */

void DiagStress::compute()
{
  int ND = fft->dimension;
  if(output->comp_strain == 0){
    for(int i=0;i<ND;i++){
      for(int j=0;j<ND;j++){
        fft->local_sigma[i][j] = 0.0;
        fft->local_sigma_NoExt[i][j] = 0.0;
        fft->avepst[i][j] = 0.0;
        fft->avepsts[i][j] = 0.0;
        fft->avepsd[i][j] = 0.0;
        fft->aveps[i][j] = 0.0;
        fft->ave_eps[i][j] = 0.0;
        fft->ave_sigma[i][j] = 0.0;
        fft->aveps_NoExt[i][j] = 0.0;
        fft->ave_eps_NoExt[i][j] = 0.0;
        fft->ave_sigma_NoExt[i][j] = 0.0;
        if (app->nphase>1) {fft->aveVs[i][j] = 0.0;}
      }
    }
    fft->strain();
    fft->total_average_strain();
    fft->stress();
    output->comp_strain = 1;
  }
  else if(output->comp_stress == 0){
    for(int i=0;i<ND;i++){
      for(int j=0;j<ND;j++){
        fft->local_sigma[i][j] = 0.0;
        fft->local_sigma_NoExt[i][j] = 0.0;
        fft->ave_sigma[i][j] = 0.0;
        fft->ave_sigma_NoExt[i][j] = 0.0;
      }
    }
    fft->stress();
    output->comp_stress = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DiagStress::stats(char *strtmp)
{
  sprintf(strtmp," %8g %8g %8g %8g %8g %8g %8g %8g %8g %8g %8g %8g",fft->ave_sigma[0][0],fft->ave_sigma[1][1],fft->ave_sigma[2][2],fft->ave_sigma[0][1],fft->ave_sigma[0][2],fft->ave_sigma[1][2],
  fft->ave_sigma_NoExt[0][0],fft->ave_sigma_NoExt[1][1],fft->ave_sigma_NoExt[2][2],fft->ave_sigma_NoExt[0][1],fft->ave_sigma_NoExt[0][2],fft->ave_sigma_NoExt[1][2]);
}

/* ---------------------------------------------------------------------- */

void DiagStress::stats_header(char *strtmp)
{
  sprintf(strtmp," %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s","s11","s22","s33","s12","s13","s23","s11NE","s22NE","s33NE","s12NE","s13NE","s23NE");
}
