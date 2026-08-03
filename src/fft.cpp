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
  condnorder = 1;
  conductivity = 0;
  xi = NULL;
  xi_sum = NULL;
  
  condxi = NULL;
  xinorm = 0;
  errorconv = 100;
  convdot = 0;
  dreg = NULL;
  dreg_sum = NULL;
  vstrain = NULL;
  phase = NULL;
  xflag = 0;
  peierlsz = 0;
  centerz1 = 0;
  centerz2 = 0;
  centerx1 = 0;
  centerx2 = 0;
  peierlsdreg = 0.0;
  X1flag = 0.0;
  X2flag = 0.0;
  Hflag = 0.0;
  eff_vol = nx*ny*nz;
  intave = 0;
  latave = 0;
  extave = 0;

  fx = fy = fz = NULL;
  f = r = NULL;
  AA = BB = FF = PD = PD2 = PDHet = DD = EE = NULL;
  MM = NULL;
  // heterogeneous form
  AAhet = BBhet = DDhet = FFhet = EEhet = NULL;
  AAhet2 = BBhet2 = DDhet2 = FFhet2 = EEhet2 = NULL;
  MGreen = NULL;
  //Grad
  gradx = grady = gradz = NULL;
  theta = NULL;
  hess = NULL;
  nhn = NULL;
  picked = NULL;

  fcore = df1core = df2core = df3core = dE_core = dE_gb = NULL;
  data_sigma = data_asigma = data_tsigma = data_psigma = data_usigma = data_rsigma = tau = NULL;

  //Grad FCC
  gradx_d = grady_d = gradz_d = NULL;
  theta_d = NULL;
  hess_d = NULL;
  nhn_d = NULL;
  picked_d = NULL;
  dt_kmc = NULL;

  datareal = NULL;
  E_int = NULL;
  E_ext = NULL;
  E_lat = NULL;
  sigma_rot = NULL;
  data_sigma_NoExt = NULL;
  xn = xb = NULL;
  sigma = NULL;
  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
  delta = ddelta = NULL;

  Rij = NULL;
  Fij = NULL;
  Fijinv = NULL;

  Gamma = NULL;

  box_exist = 0;

  lattice = NULL;
  //~ for (int i=0;i<5;i++) {material[i] = NULL;};
  material = NULL;
  materialA = NULL;
  materialB = NULL;
  materialC = NULL;
  materialD = NULL;
  materialE = NULL;
  material2 = NULL;
  material3 = NULL;
  dislocation = NULL;
  nregion = maxregion = 0;
  regions = NULL;

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      KCK[i][j] = 0.;
      KCK_CTB[i][j] = 0.;
    }
  }

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

  ortrot[0][0] = 1;
  ortrot[1][0] = 0;
  ortrot[2][0] = 0;
  ortrot[0][1] = 0;
  ortrot[1][1] = 1;
  ortrot[2][1] = 0;
  ortrot[0][2] = 0;
  ortrot[1][2] = 0;
  ortrot[2][2] = 1;

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
      // Defining the orthogonal matrix using Gram-Schmidt orthogonalization
      // defining e_1 (already normalized)
      // double inner1, inner2;
      // ortrot[0][0] = prim[0][0];
      // ortrot[1][0] = prim[1][0];
      // ortrot[2][0] = prim[2][0];

      // ortrot[0][0] = 1.0;
      // ortrot[1][0] = -1.0;
      // ortrot[2][0] = 0.0;
      // inner1 = sqrt(ortrot[0][0]*ortrot[0][0] + ortrot[1][0]*ortrot[1][0] + ortrot[2][0]*ortrot[2][0]);
      // ortrot[0][0] /= inner1;
      // ortrot[1][0] /= inner1;
      // ortrot[2][0] /= inner1;

      // // u_2 = v_2 - proj_u1 (v_2) ; e_2 = u_2/||u_2||
      // inner1 = ortrot[0][0]*prim[0][1] + ortrot[1][0]*prim[1][1] + ortrot[2][0]*prim[2][1];
      // inner2 = ortrot[0][0]*ortrot[0][0] + ortrot[1][0]*ortrot[1][0] + ortrot[2][0]*ortrot[2][0];
      // ortrot[0][1] = prim[0][1] - inner1*ortrot[0][0]/inner2;
      // ortrot[1][1] = prim[1][1] - inner1*ortrot[1][0]/inner2;
      // ortrot[2][1] = prim[2][1] - inner1*ortrot[2][0]/inner2;
      // // test: rotation using the dislocation coordinate system
      // ortrot[0][1] = 1.0;
      // ortrot[1][1] = 1.0;
      // ortrot[2][1] = -2.0;
      // // normalize e_2 = u_2/||u_2||
      // inner1 = sqrt(ortrot[0][1]*ortrot[0][1] + ortrot[1][1]*ortrot[1][1] + ortrot[2][1]*ortrot[2][1]);
      // ortrot[0][1] /= inner1;
      // ortrot[1][1] /= inner1;
      // ortrot[2][1] /= inner1;

      // // u_3 = v_3 - proj_u1 (v_3) - proj_u2 (v_3) ; e_3 = u_3/||u_3||
      // inner1 = ortrot[0][0]*prim[0][2] + ortrot[1][0]*prim[1][2] + ortrot[2][0]*prim[2][2];
      // inner2 = ortrot[0][0]*ortrot[0][0] + ortrot[1][0]*ortrot[1][0] + ortrot[2][0]*ortrot[2][0];
      // ortrot[0][2] = prim[0][2] - inner1*ortrot[0][0]/inner2;
      // ortrot[1][2] = prim[1][2] - inner1*ortrot[1][0]/inner2;
      // ortrot[2][2] = prim[2][2] - inner1*ortrot[2][0]/inner2;
      // inner1 = ortrot[0][1]*prim[0][2] + ortrot[1][1]*prim[1][2] + ortrot[2][1]*prim[2][2];
      // inner2 = ortrot[0][1]*ortrot[0][1] + ortrot[1][1]*ortrot[1][1] + ortrot[2][1]*ortrot[2][1];
      // ortrot[0][2] -= inner1*ortrot[0][1]/inner2;
      // ortrot[1][2] -= inner1*ortrot[1][1]/inner2;
      // ortrot[2][2] -= inner1*ortrot[2][1]/inner2;
      // // test: rotation using the dislocation coordinate system
      // ortrot[0][2] = 1.0;
      // ortrot[1][2] = 1.0;
      // ortrot[2][2] = 1.0;
      // // normalize e_2 = u_2/||u_2||
      // inner1 = sqrt(ortrot[0][2]*ortrot[0][2] + ortrot[1][2]*ortrot[1][2] + ortrot[2][2]*ortrot[2][2]);
      // ortrot[0][2] /= inner1;
      // ortrot[1][2] /= inner1;
      // ortrot[2][2] /= inner1;
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

  // Calculate the inverse of the ortrot
  norvolort = ortrot[1][0]*ortrot[2][1]*ortrot[0][2] - ortrot[2][0]*ortrot[1][1]*ortrot[0][2] +
        ortrot[2][0]*ortrot[0][1]*ortrot[1][2] - ortrot[0][0]*ortrot[2][1]*ortrot[1][2] +
        ortrot[0][0]*ortrot[1][1]*ortrot[2][2] - ortrot[1][0]*ortrot[0][1]*ortrot[2][2];

  invortrot[0][0] = (ortrot[1][1]*ortrot[2][2] - ortrot[2][1]*ortrot[1][2])/norvolort;
  invortrot[0][1] = (ortrot[2][1]*ortrot[0][2] - ortrot[0][1]*ortrot[2][2])/norvolort;
  invortrot[0][2] = (ortrot[0][1]*ortrot[1][2] - ortrot[1][1]*ortrot[0][2])/norvolort;

  invortrot[1][0] = (ortrot[1][2]*ortrot[2][0] - ortrot[2][2]*ortrot[1][0])/norvolort;
  invortrot[1][1] = (ortrot[2][2]*ortrot[0][0] - ortrot[0][2]*ortrot[2][0])/norvolort;
  invortrot[1][2] = (ortrot[0][2]*ortrot[1][0] - ortrot[1][2]*ortrot[0][0])/norvolort;

  invortrot[2][0] = (ortrot[1][0]*ortrot[2][1] - ortrot[2][0]*ortrot[1][1])/norvolort;
  invortrot[2][1] = (ortrot[2][0]*ortrot[0][1] - ortrot[0][0]*ortrot[2][1])/norvolort;
  invortrot[2][2] = (ortrot[0][0]*ortrot[1][1] - ortrot[1][0]*ortrot[0][1])/norvolort;
  
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
  //~ if (material) delete material;
  //~ material = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
  if (Nmaterials == 0) {
    materialA = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
    material = materialA;
  }
  else if (Nmaterials == 1) {
    materialB = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
    //~ material = materialB;
  }
  else if (Nmaterials == 2) {
    materialC = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
    //~ material = materialC;
  }
  else if (Nmaterials == 3) {
    materialD = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
    //~ material = materialD;
  }
  else if (Nmaterials ==4 ) {
    materialE = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
    //~ material = materialE;
  }
  if (Nmaterials < 5) {
    Nmaterials += 1; // increase the counter 
  }
}

// New material and region input
void FFT::set_material2(int narg, char **arg)
{
  if (material2) delete material2;
  material2 = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
}

void FFT::set_material3(int narg, char **arg)
{
  if (material3) delete material3;
  material3 = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
}

void FFT::set_dislocation(int narg, char **arg)
{
  if (dislocation) delete dislocation;
  dislocation = new Material(pfdd_p,narg,arg,nx,ny,nz,app->slip_systems,app->num_planes);
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
/* ----------------------------------------------------------------------
   Calculate determinant o 6x6 tensor
------------------------------------------------------------------------- */
double FFT::determinant(double a[6][6], int k)
{
  double s, det;
    double b[6][6];
    int i, j, m, n, c;

    s = 1.;
    det = 0.;

    if (k == 1) {
      return (a[0][0]);
    }
    else {
      det = 0;
      for (c = 0; c < k; c++) {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++) {
          for (j = 0 ;j < k; j++) {
            b[i][j] = 0.;
            if (i != 0 && j != c){
              b[m][n] = a[i][j];
              if (n < (k - 2)) {
                n++;
              }
              else {
                n = 0;
                m++;
              }
            }
          }
        }
        det = det + s * (a[0][c] * determinant(b, k - 1));
        s = -1. * s;
      }
    }
    return (det);
}

void FFT::rotate_stiffness(double CI[6][6], double ed1[3], double ed2[3], double ed3[3])
{
  double omega[6][6], K[6][6], KT[6][6], KC[6][6]; //, KCK[6][6];
  double sum_temp;
  int i, j, k;
  double e1[3], e2[3], e3[3];

  for (i = 0; i < 3; i++) {
    e1[i] = 0.;
    e2[i] = 0.;
    e3[i] = 0.;
  }
  e1[0] = 1.0; e2[1] = 1.0; e3[2] = 1.0;

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      omega[i][j] = 0.;
    }
  }
  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      K[i][j] = 0.;
      KT[i][j] = 0.;
      KC[i][j] = 0.;
      // KCK[i][j] = 0.;
    }
  }

  // omega

  // omega[0][0] = ed1[0]*e1[0] + ed1[1]*e1[1] + ed1[2]*e1[2];
  // omega[0][1] = ed2[0]*e1[0] + ed2[1]*e1[1] + ed2[2]*e1[2];
  // omega[0][2] = ed3[0]*e1[0] + ed3[1]*e1[1] + ed3[2]*e1[2];
  //
  // omega[1][0] = ed1[0]*e2[0] + ed1[1]*e2[1] + ed1[2]*e2[2];
  // omega[1][1] = ed2[0]*e2[0] + ed2[1]*e2[1] + ed2[2]*e2[2];
  // omega[1][2] = ed3[0]*e2[0] + ed3[1]*e2[1] + ed3[2]*e2[2];
  //
  // omega[2][0] = ed1[0]*e3[0] + ed1[1]*e3[1] + ed1[2]*e3[2];
  // omega[2][1] = ed2[0]*e3[0] + ed2[1]*e3[1] + ed2[2]*e3[2];
  // omega[2][2] = ed3[0]*e3[0] + ed3[1]*e3[1] + ed3[2]*e3[2];

  omega[0][0] = ed1[0];
  omega[0][1] = ed1[1];
  omega[0][2] = ed1[2];

  omega[1][0] = ed2[0];
  omega[1][1] = ed2[1];
  omega[1][2] = ed2[2];

  omega[2][0] = ed3[0];
  omega[2][1] = ed3[1];
  omega[2][2] = ed3[2];

  // K1

  K[0][0] = pow(omega[0][0], 2.);
  K[1][0] = pow(omega[1][0], 2.);
  K[2][0] = pow(omega[2][0], 2.);

  K[0][1] = pow(omega[0][1], 2.);
  K[1][1] = pow(omega[1][1], 2.);
  K[2][1] = pow(omega[2][1], 2.);

  K[0][2] = pow(omega[0][2], 2.);
  K[1][2] = pow(omega[1][2], 2.);
  K[2][2] = pow(omega[2][2], 2.);

  // K2

  K[0][3] = 2. * omega[0][1] * omega[0][2];
  K[1][3] = 2. * omega[1][1] * omega[1][2];
  K[2][3] = 2. * omega[2][1] * omega[2][2];

  K[0][4] = 2. * omega[0][2] * omega[0][0];
  K[1][4] = 2. * omega[1][2] * omega[1][0];
  K[2][4] = 2. * omega[2][2] * omega[2][0];

  K[0][5] = 2. * omega[0][0] * omega[0][1];
  K[1][5] = 2. * omega[1][0] * omega[1][1];
  K[2][5] = 2. * omega[2][0] * omega[2][1];

  // K3

  K[3][0] = omega[1][0] * omega[2][0];
  K[4][0] = omega[2][0] * omega[0][0];
  K[5][0] = omega[0][0] * omega[1][0];

  K[3][1] = omega[1][1] * omega[2][1];
  K[4][1] = omega[2][1] * omega[0][1];
  K[5][1] = omega[0][1] * omega[1][1];

  K[3][2] = omega[1][2] * omega[2][2];
  K[4][2] = omega[2][2] * omega[0][2];
  K[5][2] = omega[0][2] * omega[1][2];

  // K4

  K[3][3] = omega[1][1] * omega[2][2] + omega[1][2] * omega[2][1];
  K[4][3] = omega[2][1] * omega[0][2] + omega[2][2] * omega[0][1];
  K[5][3] = omega[0][1] * omega[1][2] + omega[0][2] * omega[1][1];

  K[3][4] = omega[1][2] * omega[2][0] + omega[1][0] * omega[2][2];
  K[4][4] = omega[2][2] * omega[0][0] + omega[2][0] * omega[0][2];
  K[5][4] = omega[0][2] * omega[1][0] + omega[0][0] * omega[1][2];

  K[3][5] = omega[1][0] * omega[2][1] + omega[1][1] * omega[2][0];
  K[4][5] = omega[2][0] * omega[0][1] + omega[2][1] * omega[0][0];
  K[5][5] = omega[0][0] * omega[1][1] + omega[0][1] * omega[1][0];

  // KT: transpose of K

  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      KT[j][i] = K[i][j];
    }
  }
  // KC = K * CI

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {

      sum_temp = 0.;
      for (k = 0; k < 6; k++) {
        sum_temp = sum_temp + K[i][k] * CI[k][j];
      }
      KC[i][j] = sum_temp;

    }
  }

  // KCK = KC * KT

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {

      sum_temp = 0.;
      for (k = 0; k < 6; k++) {
        sum_temp = sum_temp + KC[i][k] * KT[k][j];
      }
      KCK[i][j] = sum_temp;
      if(me == 0){
        printf("KCK[%d][%d] = %E\n", i, j, KCK[i][j]);
      }

    }
  }
}

void FFT::rotate_stiffness_CTB(double CI[6][6], double ed1[3], double ed2[3], double ed3[3])
{
  double omega[6][6], K[6][6], KT[6][6], KC[6][6]; //, KCK[6][6];
  double sum_temp;
  int i, j, k;
  double e1[3], e2[3], e3[3];

  for (i = 0; i < 3; i++) {
    e1[i] = 0.;
    e2[i] = 0.;
    e3[i] = 0.;
  }
  e1[0] = 1.0; e2[1] = 1.0; e3[2] = 1.0;

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      omega[i][j] = 0.;
    }
  }
  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      K[i][j] = 0.;
      KT[i][j] = 0.;
      KC[i][j] = 0.;
      // KCK[i][j] = 0.;
    }
  }

  // omega

  // omega[0][0] = ed1[0]*e1[0] + ed1[1]*e1[1] + ed1[2]*e1[2];
  // omega[0][1] = ed2[0]*e1[0] + ed2[1]*e1[1] + ed2[2]*e1[2];
  // omega[0][2] = ed3[0]*e1[0] + ed3[1]*e1[1] + ed3[2]*e1[2];
  //
  // omega[1][0] = ed1[0]*e2[0] + ed1[1]*e2[1] + ed1[2]*e2[2];
  // omega[1][1] = ed2[0]*e2[0] + ed2[1]*e2[1] + ed2[2]*e2[2];
  // omega[1][2] = ed3[0]*e2[0] + ed3[1]*e2[1] + ed3[2]*e2[2];
  //
  // omega[2][0] = ed1[0]*e3[0] + ed1[1]*e3[1] + ed1[2]*e3[2];
  // omega[2][1] = ed2[0]*e3[0] + ed2[1]*e3[1] + ed2[2]*e3[2];
  // omega[2][2] = ed3[0]*e3[0] + ed3[1]*e3[1] + ed3[2]*e3[2];

  omega[0][0] = ed1[0];
  omega[0][1] = ed1[1];
  omega[0][2] = ed1[2];

  omega[1][0] = ed2[0];
  omega[1][1] = ed2[1];
  omega[1][2] = ed2[2];

  omega[2][0] = ed3[0];
  omega[2][1] = ed3[1];
  omega[2][2] = ed3[2];

  // K1

  K[0][0] = pow(omega[0][0], 2.);
  K[1][0] = pow(omega[1][0], 2.);
  K[2][0] = pow(omega[2][0], 2.);

  K[0][1] = pow(omega[0][1], 2.);
  K[1][1] = pow(omega[1][1], 2.);
  K[2][1] = pow(omega[2][1], 2.);

  K[0][2] = pow(omega[0][2], 2.);
  K[1][2] = pow(omega[1][2], 2.);
  K[2][2] = pow(omega[2][2], 2.);

  // K2

  K[0][3] = 2. * omega[0][1] * omega[0][2];
  K[1][3] = 2. * omega[1][1] * omega[1][2];
  K[2][3] = 2. * omega[2][1] * omega[2][2];

  K[0][4] = 2. * omega[0][2] * omega[0][0];
  K[1][4] = 2. * omega[1][2] * omega[1][0];
  K[2][4] = 2. * omega[2][2] * omega[2][0];

  K[0][5] = 2. * omega[0][0] * omega[0][1];
  K[1][5] = 2. * omega[1][0] * omega[1][1];
  K[2][5] = 2. * omega[2][0] * omega[2][1];

  // K3

  K[3][0] = omega[1][0] * omega[2][0];
  K[4][0] = omega[2][0] * omega[0][0];
  K[5][0] = omega[0][0] * omega[1][0];

  K[3][1] = omega[1][1] * omega[2][1];
  K[4][1] = omega[2][1] * omega[0][1];
  K[5][1] = omega[0][1] * omega[1][1];

  K[3][2] = omega[1][2] * omega[2][2];
  K[4][2] = omega[2][2] * omega[0][2];
  K[5][2] = omega[0][2] * omega[1][2];

  // K4

  K[3][3] = omega[1][1] * omega[2][2] + omega[1][2] * omega[2][1];
  K[4][3] = omega[2][1] * omega[0][2] + omega[2][2] * omega[0][1];
  K[5][3] = omega[0][1] * omega[1][2] + omega[0][2] * omega[1][1];

  K[3][4] = omega[1][2] * omega[2][0] + omega[1][0] * omega[2][2];
  K[4][4] = omega[2][2] * omega[0][0] + omega[2][0] * omega[0][2];
  K[5][4] = omega[0][2] * omega[1][0] + omega[0][0] * omega[1][2];

  K[3][5] = omega[1][0] * omega[2][1] + omega[1][1] * omega[2][0];
  K[4][5] = omega[2][0] * omega[0][1] + omega[2][1] * omega[0][0];
  K[5][5] = omega[0][0] * omega[1][1] + omega[0][1] * omega[1][0];

  // KT: transpose of K

  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      KT[j][i] = K[i][j];
    }
  }
  // KC = K * CI

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {

      sum_temp = 0.;
      for (k = 0; k < 6; k++) {
        sum_temp = sum_temp + K[i][k] * CI[k][j];
      }
      KC[i][j] = sum_temp;

    }
  }

  // KCK = KC * KT

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {

      sum_temp = 0.;
      for (k = 0; k < 6; k++) {
        sum_temp = sum_temp + KC[i][k] * KT[k][j];
      }
      KCK_CTB[i][j] = sum_temp;
      if(me == 0){
        printf("KCK_CTB[%d][%d] = %E\n", i, j, KCK_CTB[i][j]);
      }

    }
  }
}
void FFT::rotate_stiffnessfour(double Cfour[3][3][3][3], double R[3][3])
{
  // rotate 4th order stiffness tensor
  int i,j,k,l,m,n,o,p;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {
          RCR[i][j][k][l] = 0.0;
        }
      }
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {
          for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
              for (o = 0; o < 3; o++) {
                for (p = 0; p < 3; p++) {
                  // r_mi r_nj r_ok r_ pl C_ijkl
                  RCR[m][n][o][p] += R[m][i]*R[n][j]*R[o][k]*R[p][l]*Cfour[i][j][k][l];
                }
              }
            }
          }
        }
      }
    }
  }

  if(me == 0){
    printf("--------------");
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          for (l = 0; l < 3; l++) {
            printf("rotated = %E\n", RCR[i][j][k][l]);
          }
        }
      }
    }
  }

}
