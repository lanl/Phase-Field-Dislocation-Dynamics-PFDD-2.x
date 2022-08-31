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

#include "math.h"
#include "mpi.h"
#include <stdio.h>
#include "stdlib.h"
#include "string.h"
#include "fft.h"
#include "app.h"
#include "lattice.h"
#include "material.h"
#include "memory.h"
#include "error.h"

#include "style_region.h"

using namespace PFDD_NS;

#define DELTA 1

/* ---------------------------------------------------------------------- */

FFT::FFT(PFDD_C *pfdd_p, int narg, char **arg) : Pointers(pfdd_p)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  dimension = 3;
  xperiodic = yperiodic = zperiodic = 1;
  nonperiodic = 0;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nx = ny = nz = 0;
  norder = 1;
  xi = NULL;
  xi_sum = NULL;

  fx = fy = fz = NULL;
  f = r = NULL;
  BB = FF = DD = NULL;
  //Grad
  gradx = grady = gradz = NULL;
  theta = NULL;

  fcore = df1core = df2core = df3core = dE_core = NULL;
  data_sigma = tau = NULL;
  xn = xb = NULL;
  sigma = NULL;
  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
  delta = ddelta = NULL;

  box_exist = 0;

  lattice = NULL;
  material = NULL;
  nregion = maxregion = 0;
  regions = NULL;

  primitive = 0;

  prim[0][0] = 1;
  prim[1][0] = 0;
  prim[2][0] = 0;
  prim[0][1] = 0;
  prim[1][1] = 1;
  prim[2][1] = 0;
  prim[0][2] = 0;
  prim[1][2] = 0;
  prim[2][2] = 1;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"primitive") == 0) {
      if (iarg+12 > narg) error->all(FLERR,"Illegal fft command");
      primitive = 1;
      iarg += 1;
      if (strcmp(arg[iarg],"x") == 0) {
	iarg += 1;
	prim[0][0] = atof(arg[iarg]);
	prim[1][0] = atof(arg[iarg+1]);
	prim[2][0] = atof(arg[iarg+2]);
	double norm = sqrt(prim[0][0]*prim[0][0]+prim[1][0]*prim[1][0]+prim[2][0]*prim[2][0]);
	prim[0][0] /= norm;
	prim[1][0] /= norm;
	prim[2][0] /= norm;
	iarg += 3;
      }
      if (strcmp(arg[iarg],"y") == 0) {
	iarg += 1;
	prim[0][1] = atof(arg[iarg]);
	prim[1][1] = atof(arg[iarg+1]);
	prim[2][1] = atof(arg[iarg+2]);
	double norm = sqrt(prim[0][1]*prim[0][1]+prim[1][1]*prim[1][1]+prim[2][1]*prim[2][1]);
	prim[0][1] /= norm;
	prim[1][1] /= norm;
	prim[2][1] /= norm;
	iarg += 3;
      }
      if (strcmp(arg[iarg],"z") == 0) {
	iarg += 1;
	prim[0][2] = atof(arg[iarg]);
	prim[1][2] = atof(arg[iarg+1]);
	prim[2][2] = atof(arg[iarg+2]);
	double norm = sqrt(prim[0][2]*prim[0][2]+prim[1][2]*prim[1][2]+prim[2][2]*prim[2][2]);
	prim[0][2] /= norm;
	prim[1][2] /= norm;
	prim[2][2] /= norm;
	iarg += 3;
      }
    }
    else if (strcmp(arg[iarg],"pbc") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fft command");
      iarg += 1;
      xperiodic = atoi(arg[iarg]);
      yperiodic = atoi(arg[iarg+1]);
      zperiodic = atoi(arg[iarg+2]);
      periodicity[0] = xperiodic;
      periodicity[1] = yperiodic;
      periodicity[2] = zperiodic;
      iarg += 3;
    }
    else if (strcmp(arg[iarg],"mode") == 0) {
      iarg += 1;
      mode = atoi(arg[iarg]);
      iarg += 1;
    }
    else error->all(FLERR,"Illegal fft_style command");
  }

  // Calculate the inverse of the prim
  norvol = prim[1][0]*prim[2][1]*prim[0][2] - prim[2][0]*prim[1][1]*prim[0][2] +
        prim[2][0]*prim[0][1]*prim[1][2] - prim[0][0]*prim[2][1]*prim[1][2] +
        prim[0][0]*prim[1][1]*prim[2][2] - prim[1][0]*prim[0][1]*prim[2][2];

  invprim[0][0] = (prim[1][1]*prim[2][2] - prim[2][1]*prim[1][2])/norvol;
  invprim[0][1] = (prim[2][1]*prim[0][2] - prim[0][1]*prim[2][2])/norvol;
  invprim[0][2] = (prim[0][1]*prim[1][2] - prim[1][1]*prim[0][2])/norvol;

  invprim[1][0] = (prim[1][2]*prim[2][0] - prim[2][2]*prim[1][0])/norvol;
  invprim[1][1] = (prim[2][2]*prim[0][0] - prim[0][2]*prim[2][0])/norvol;
  invprim[1][2] = (prim[0][2]*prim[1][0] - prim[1][2]*prim[0][0])/norvol;

  invprim[2][0] = (prim[1][0]*prim[2][1] - prim[2][0]*prim[1][1])/norvol;
  invprim[2][1] = (prim[2][0]*prim[0][1] - prim[0][0]*prim[2][1])/norvol;
  invprim[2][2] = (prim[0][0]*prim[1][1] - prim[1][0]*prim[0][1])/norvol;

}

/* ---------------------------------------------------------------------- */

FFT::~FFT()
{
  delete lattice;
  for (int i = 0; i < nregion; i++) delete regions[i];
  memory->sfree(regions);
}

/* ----------------------------------------------------------------------
   setup global box
   assumes boxlo/hi are already set
------------------------------------------------------------------------- */

void FFT::set_box()
{
  if (boxxlo >= boxxhi || boxylo >= boxyhi || boxzlo >= boxzhi)
    error->one(FLERR,"Box bounds are invalid");

  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;
}



/* ----------------------------------------------------------------------
   create a lattice
   delete it if style = none
------------------------------------------------------------------------- */

void FFT::set_lattice(int narg, char **arg)
{
  if (lattice) delete lattice;
  lattice = new Lattice(pfdd_p,narg,arg);
  if (lattice->style == 0) {
    delete lattice;
    lattice = NULL;
  }
}

/* ----------------------------------------------------------------------
   create a material
   delete it if style = none
------------------------------------------------------------------------- */

void FFT::set_material(int narg, char **arg)
{
  if (material) delete material;
  material = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems);
}

/* ----------------------------------------------------------------------
   create a new region
------------------------------------------------------------------------- */

void FFT::add_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal region command");

  if (find_region(arg[0]) >= 0) error->all(FLERR,"Reuse of region ID");

  // extend Region list if necessary

  if (nregion == maxregion) {
    maxregion += DELTA;
    regions = (Region **)
      memory->srealloc(regions,maxregion*sizeof(Region *),"domain:regions");
  }

  // create the Region

  if (strcmp(arg[1],"none") == 0) error->all(FLERR,"Invalid region style");

#define REGION_CLASS
#define RegionStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) \
    regions[nregion] = new Class(pfdd_p,narg,arg);
#include "style_region.h"
#undef REGION_CLASS

  else error->all(FLERR,"Invalid region style");

  nregion++;
}

/* ----------------------------------------------------------------------
   return region index if name matches existing region ID
   return -1 if no such region
------------------------------------------------------------------------- */

int FFT::find_region(char *name)
{
  for (int iregion = 0; iregion < nregion; iregion++)
    if (strcmp(name,regions[iregion]->id) == 0) return iregion;
  return -1;
}

/* ----------------------------------------------------------------------
   boundary settings from the input script
------------------------------------------------------------------------- */

void FFT::set_boundary(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal boundary command");

  if (strcmp(arg[0],"n") == 0) xperiodic = 0;
  else if (strcmp(arg[0],"p") == 0) xperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");
  if (strcmp(arg[1],"n") == 0) yperiodic = 0;
  else if (strcmp(arg[1],"p") == 0) yperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");
  if (strcmp(arg[2],"n") == 0) zperiodic = 0;
  else if (strcmp(arg[2],"p") == 0) zperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");

  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) nonperiodic = 1;

  if (nonperiodic && app->appclass != App::LATTICE)
    error->all(FLERR,"Boundary command currently only supported by on-lattice apps");
}

/* ----------------------------------------------------------------------
   assign nprocs to 1d box as equal partitions
------------------------------------------------------------------------- */

void FFT::procs2domain_1d()
{
  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    if (user_procgrid[1] != 1 || user_procgrid[2] != 1)
      error->all(FLERR,"App style proc count is not valid for 1d simulation");
    procgrid[0] = user_procgrid[0];
  } else {
    procgrid[0] = nprocs;
  }

  procgrid[1] = procgrid[2] = 1;

  myloc[0] = me;
  myloc[1] = myloc[2] = 0;

  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1)
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo;
  subyhi = boxyhi;
  subzlo = boxzlo;
  subzhi = boxzhi;
}

/* ----------------------------------------------------------------------
   assign nprocs to 2d box so as to minimize perimeter per proc
------------------------------------------------------------------------- */

void FFT::procs2domain_2d()
{
  int ipx,ipy;
  double boxx,boxy,surf;

  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    if (user_procgrid[2] != 1)
      error->all(FLERR,"App style proc count is not valid for 2d simulation");
    procgrid[0] = user_procgrid[0];
    procgrid[1] = user_procgrid[1];

  } else {

    // loop thru all possible factorizations of nprocs
    // surf = perimeter of a proc sub-domain

    double bestsurf = 2.0 * (xprd+yprd);

    ipx = 1;
    while (ipx <= nprocs) {
      if (nprocs % ipx == 0) {
	ipy = nprocs/ipx;
	boxx = xprd/ipx;
	boxy = yprd/ipy;
	surf = boxx + boxy;
	if (surf < bestsurf) {
	  bestsurf = surf;
	  procgrid[0] = ipx;
	  procgrid[1] = ipy;
	}
      }
      ipx++;
    }
  }

  procgrid[2] = 1;

  myloc[0] = me % procgrid[0];
  myloc[1] = me/procgrid[0];
  myloc[2] = 0;

  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1)
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo + myloc[1] * yprd/procgrid[1];
  if (myloc[1] < procgrid[1]-1)
    subyhi = boxylo + (myloc[1]+1) * yprd/procgrid[1];
  else subyhi = boxyhi;

  subzlo = boxzlo;
  subzhi = boxzhi;
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d box so as to minimize surface area per proc
------------------------------------------------------------------------- */

void FFT::procs2domain_3d()
{
  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;

  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    procgrid[0] = user_procgrid[0];
    procgrid[1] = user_procgrid[1];
    procgrid[2] = user_procgrid[2];

  } else {

    double bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd);

    // loop thru all possible factorizations of nprocs
    // surf = surface area of a proc sub-domain

    ipx = 1;
    while (ipx <= nprocs) {
      if (nprocs % ipx == 0) {
	nremain = nprocs/ipx;
	ipy = 1;
	while (ipy <= nremain) {
	  if (nremain % ipy == 0) {
	    ipz = nremain/ipy;
	    boxx = xprd/ipx;
	    boxy = yprd/ipy;
	    boxz = zprd/ipz;
	    surf = boxx*boxy + boxy*boxz + boxz*boxx;
	    if (surf < bestsurf) {
	      bestsurf = surf;
	      procgrid[0] = ipx;
	      procgrid[1] = ipy;
	      procgrid[2] = ipz;
	    }
	  }
	  ipy++;
	}
      }
      ipx++;
    }
  }

  myloc[0] = me % procgrid[0];
  myloc[1] = (me/procgrid[0]) % procgrid[1];
  myloc[2] = me / (procgrid[0]*procgrid[1]);

  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1)
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo + myloc[1] * yprd/procgrid[1];
  if (myloc[1] < procgrid[1]-1)
    subyhi = boxylo + (myloc[1]+1) * yprd/procgrid[1];
  else subyhi = boxyhi;

  subzlo = boxzlo + myloc[2] * zprd/procgrid[2];
  if (myloc[2] < procgrid[2]-1)
    subzhi = boxzlo + (myloc[2]+1) * zprd/procgrid[2];
  else subzhi = boxzhi;
}
