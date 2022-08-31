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
#include "region_block.h"
#include "fft.h"
#include "error.h"

using namespace PFDD_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegBlock::RegBlock(PFDD_C *pfdd_p, int narg, char **arg) : Region(pfdd_p, narg, arg)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"INF") == 0 || strcmp(arg[2],"EDGE") == 0) {
    if (fft->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[2],"INF") == 0) xlo = -BIG;
    else xlo = fft->boxxlo;
  } else xlo = xscale*atof(arg[2]);

  if (strcmp(arg[3],"INF") == 0 || strcmp(arg[3],"EDGE") == 0) {
    if (fft->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[3],"INF") == 0) xhi = BIG;
    else xhi = fft->boxxhi;
  } else xhi = xscale*atof(arg[3]);

  if (strcmp(arg[4],"INF") == 0 || strcmp(arg[4],"EDGE") == 0) {
    if (fft->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[4],"INF") == 0) ylo = -BIG;
    else ylo = fft->boxylo;
  } else ylo = yscale*atof(arg[4]);

  if (strcmp(arg[5],"INF") == 0 || strcmp(arg[5],"EDGE") == 0) {
    if (fft->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[5],"INF") == 0) yhi = BIG;
    else yhi = fft->boxyhi;
  } else yhi = yscale*atof(arg[5]);

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (fft->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
    else zlo = fft->boxzlo;
  } else zlo = zscale*atof(arg[6]);

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (fft->box_exist == 0) 
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[7],"INF") == 0) zhi = BIG;
    else zhi = fft->boxzhi;
  } else zhi = zscale*atof(arg[7]);

  // error check

  if (xlo > xhi || ylo > yhi || zlo > zhi)
    error->all(FLERR,"Illegal region block command");

  // extent of block
  
  extent_xlo = xlo;
  extent_xhi = xhi;
  extent_ylo = ylo;
  extent_yhi = yhi;
  extent_zlo = zlo;
  extent_zhi = zhi;
}

/* ---------------------------------------------------------------------- */

int RegBlock::match(double x, double y, double z)
{
  int inside;
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    inside = 1;
  else inside = 0;

  return !(inside ^ interior);         // 1 if same, 0 if different
}
