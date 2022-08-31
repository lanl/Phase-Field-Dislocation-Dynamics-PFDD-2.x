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
#include "string.h"
#include "stdlib.h"
#include "app.h"
#include "fft.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "random_mars.h"
#include "comm_lattice.h"
#include "output.h"
#include "solve.h"

using namespace PFDD_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

App::App(PFDD_C *pfdd_p, int narg, char **arg) : Pointers(pfdd_p)
{
  if (narg < 1) error->all(FLERR,"Illegal app_style command");

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  appclass = LATTICE;
  time = 0.0;
  first_run = 1;
  dimension = 3;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  temperature = 0.0;
  ninteger = ndouble = 0;
  id = NULL;
  xyz = NULL;
  sigma = NULL;
  sa = sb = sc = 1.0;

  comm = NULL;
  ranapp = NULL;
  ranstrict = NULL;

  nlocal = nghost = nmax = 0;
  owner = NULL;
  index = NULL;

  neighlayer = 1;
  maxneigh = 0;
  numneigh = NULL;
  neighbor = NULL;
  sites_exist = 0;
  non_Schmid = 0;
  sites = NULL;

  core_flag = -1;
  setup_flag = -1;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"core_energy") == 0) {
      iarg += 1;
      core_flag = atoi(arg[iarg]);
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"setup") == 0) {
      iarg += 1;
      set_initial_sxtal(arg[iarg]);
      iarg += 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

App::~App()
{
  delete [] style;

  memory->destroy(id);
  memory->destroy(xyz);

  delete comm;
  memory->destroy(owner);
  memory->destroy(index);

  memory->destroy(numneigh);
  memory->destroy(neighbor);

}

/* ---------------------------------------------------------------------- */

void App::create_arrays()
{

}

/* ---------------------------------------------------------------------- */

void App::recreate_arrays()
{

  create_arrays();
}

/* ---------------------------------------------------------------------- */

void App::run(int narg, char **arg)
{
  if (appclass != GENERAL && fft->box_exist == 0)
    error->all(FLERR,"Cannot run application until simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal run command");

  stoptime = -1;
  stopsteps = -1;
  stoptime = time + atof(arg[0]);
  if (stoptime < time) error->all(FLERR,"Illegal run command");

  // read optional args

  int uptoflag = 0;
  int preflag = 1;
  int postflag = 1;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"upto") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal run command");
      uptoflag = 1;
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"pre") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) preflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) preflag = 1;
      else error->all(FLERR,"Illegal run command");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"post") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) postflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) postflag = 1;
      else error->all(FLERR,"Illegal run command");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"steps") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run command");
      stopsteps = atoi(arg[iarg+1]);
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"time") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run command");
      stoptime = atof(arg[iarg+1]);
      iarg += 2;
    }
    else error->all(FLERR,"Illegal run command");
  }

  // adjust stoptime if upto was specified

  if (uptoflag) stoptime -= time;
  if (uptoflag && stoptime < 0.0)
    error->all(FLERR,"Run upto value is before current time");

  // perform a single run via app's init(), setup(), and iterate()
  // if pre or 1st run, do app init
  // setup computes initial propensities
  // if post, do full Finish, else just print time

  if (preflag || first_run) {
    init();
    first_run = 0;
  }

  if (fft->me == 0) {
    if (screen) fprintf(screen,"Setting up run ...\n");
    if (logfile) fprintf(logfile,"Setting up run ...\n");
  }

  timer->init();
  setup();
  if (stoptime > time) iterate();

  Finish finish(pfdd_p,postflag);
}

/* ---------------------------------------------------------------------- */

void App::reset_time(double newtime)
{
  time = newtime;
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   if don't recognize name, pass it along to lower-level app
   names iarrayN and orderparamN mean entry N from 1 to ninteger or ndouble
 ------------------------------------------------------------------------- */

void *App::extract(char *name)
{
  if (strcmp(name,"dimension") == 0) return (void *) &fft->dimension;
  if (strcmp(name,"boxxlo") == 0) return (void *) &fft->boxxlo;
  if (strcmp(name,"boxxhi") == 0) return (void *) &fft->boxxhi;
  if (strcmp(name,"boxylo") == 0) return (void *) &fft->boxylo;
  if (strcmp(name,"boxyhi") == 0) return (void *) &fft->boxyhi;
  if (strcmp(name,"boxzlo") == 0) return (void *) &fft->boxzlo;
  if (strcmp(name,"boxzhi") == 0) return (void *) &fft->boxzhi;

  if (strcmp(name,"nglobal") == 0) return (void *) &nglobal;
  if (strcmp(name,"nlocal") == 0) return (void *) &nlocal;

  if (strcmp(name,"id") == 0) return (void *) id;
  if (strcmp(name,"xyz") == 0) return (void *) xyz;

  if (strcmp(name,"site") == 0) {
    if (ninteger == 0) return NULL;
    return (void *) fft->xi[0];
  }

  // if (strstr(name,"xi") == name) {
  //   int n = atoi(&name[6]);
  //   if (n < 1 || n > ndouble) return NULL;
  //   return (void *) xi[0][n-1];
  // }

  return extract_app(name);
}

/* ----------------------------------------------------------------------
   return max ID
   may not be nglobal if site IDs are not contiguous
 ------------------------------------------------------------------------- */

tagint App::max_site_ID()
{
  tagint max = -1;
  for (int i = 0; i < nlocal; i++) max = MAX(max,id[i]);
  tagint all;
  MPI_Allreduce(&max,&all,1,MPI_PFDD_TAGINT,MPI_MAX,world);
  return all;
}

/* ----------------------------------------------------------------------
   check if site IDs are contiguous from 1 to N
 ------------------------------------------------------------------------- */

int App::contiguous_sites()
{
  tagint min = nglobal+1;
  tagint max = -1;

  for (int i = 0; i < nlocal; i++) {
    min = MIN(min,id[i]);
    max = MAX(max,id[i]);
  }

  tagint all;
  MPI_Allreduce(&min,&all,1,MPI_PFDD_TAGINT,MPI_MIN,world);
  if (all != 1) return 0;
  MPI_Allreduce(&max,&all,1,MPI_PFDD_TAGINT,MPI_MAX,world);
  if (all != nglobal) return 0;
  return 1;
}

void App::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"update_only") == 0) set_update_only(narg,arg);
  else if (strcmp(command,"sigma") == 0) set_sigma(narg,arg);
  else if (strcmp(command,"scale") == 0) set_scale(narg,arg);
  else if (strcmp(command,"non_Schmid") == 0) set_nonSchmid(narg,arg);
  else input_app(command,narg,arg);
}

void App::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ---------------------------------------------------------------------- */

void App::set_update_only(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal update_only command");
  if (strcmp(arg[0],"yes") == 0) update_only = 1;
  else if (strcmp(arg[0],"no") == 0) update_only = 0;
  else error->all(FLERR,"Illegal update_only command");

  if (update_only && !allow_update)
    error->all(FLERR,"App does not permit user_update yes");
}

/* ----------------------------------------------------------------------
   setup dimension
------------------------------------------------------------------------- */

void App::set_dimension(int narg, char **arg)
{
  dimension = atoi(arg[0]);
}

void App::set_sigma(int narg, char **arg)
{
  if (narg <= 6) error->all(FLERR,"Illegal sigma command");

  memory->create(sigma,dimension,dimension,"sigma");
  memory->create(deltasig,dimension,dimension,"deltasig");
  for (int i=0; i<dimension; i++){
    for (int j=0; j<dimension; j++){
      sigma[i][j] = 0;
      deltasig[i][j] = 0;
    }
  }

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"initial") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal sigma command");
      iarg += 1;
      sigma[0][0] = atof(arg[iarg]);
      iarg += 1;
      sigma[1][1] = atof(arg[iarg]);
      iarg += 1;
      sigma[2][2] = atof(arg[iarg]);
      iarg += 1;
      sigma[0][1] = sigma[1][0] = atof(arg[iarg]);
      iarg += 1;
      sigma[0][2] = sigma[2][0] = atof(arg[iarg]);
      iarg += 1;
      sigma[1][2] = sigma[2][1] = atof(arg[iarg]);
      iarg += 1;
  /* for (int i=0; i<dimension; i++) {
    for (int j=0; j<dimension; j++) {
      if(me == 0){
        if (logfile) fprintf(logfile,"app->sigma[%d][%d] = %lf\n", i, j, sigma[i][j]);
        if (screen) fprintf(screen,"app->sigma[%d][%d] = %lf\n", i, j, sigma[i][j]);
      }
    }
  } */
    }
    else if (strcmp(arg[iarg],"delta") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal sigma command");
      iarg += 1;
      deltasig[0][0] = atof(arg[iarg]);
      iarg += 1;
      deltasig[1][1] = atof(arg[iarg]);
      iarg += 1;
      deltasig[2][2] = atof(arg[iarg]);
      iarg += 1;
      deltasig[0][1] = deltasig[1][0] = atof(arg[iarg]);
      iarg += 1;
      deltasig[0][2] = deltasig[2][0] = atof(arg[iarg]);
      iarg += 1;
      deltasig[1][2] = deltasig[2][1] = atof(arg[iarg]);
      iarg += 1;
    }
  }
}

void App::set_scale(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal scale command");
  sa = atof(arg[0]);
  sb = atof(arg[1]);
  sc = atof(arg[2]);

  if(me == 0){
    if (logfile) fprintf(logfile,"sa = %lf, sb = %lf, sc = %lf\n", sa, sb, sc);
    if (screen) fprintf(screen,"sa = %lf, sb = %lf, sc = %lf\n", sa, sb, sc);
  }


  fft->sclprim[0][0] = sa*fft->prim[0][0];
  fft->sclprim[1][0] = sa*fft->prim[1][0];
  fft->sclprim[2][0] = sa*fft->prim[2][0];
  fft->sclprim[0][1] = sb*fft->prim[0][1];
  fft->sclprim[1][1] = sb*fft->prim[1][1];
  fft->sclprim[2][1] = sb*fft->prim[2][1];
  fft->sclprim[0][2] = sc*fft->prim[0][2];
  fft->sclprim[1][2] = sc*fft->prim[1][2];
  fft->sclprim[2][2] = sc*fft->prim[2][2];

  // Calculate volume of lattce cell
  fft->vol = fft->sclprim[1][0]*fft->sclprim[2][1]*fft->sclprim[0][2] - fft->sclprim[2][0]*fft->sclprim[1][1]*fft->sclprim[0][2] +
             fft->sclprim[2][0]*fft->sclprim[0][1]*fft->sclprim[1][2] - fft->sclprim[0][0]*fft->sclprim[2][1]*fft->sclprim[1][2] +
             fft->sclprim[0][0]*fft->sclprim[1][1]*fft->sclprim[2][2] - fft->sclprim[1][0]*fft->sclprim[0][1]*fft->sclprim[2][2];

  fft->four[0][0] = (fft->sclprim[1][1]*fft->sclprim[2][2] - fft->sclprim[2][1]*fft->sclprim[1][2])/fft->vol;
  fft->four[0][1] = (fft->sclprim[2][1]*fft->sclprim[0][2] - fft->sclprim[0][1]*fft->sclprim[2][2])/fft->vol;
  fft->four[0][2] = (fft->sclprim[0][1]*fft->sclprim[1][2] - fft->sclprim[1][1]*fft->sclprim[0][2])/fft->vol;

  fft->four[1][0] = (fft->sclprim[1][2]*fft->sclprim[2][0] - fft->sclprim[2][2]*fft->sclprim[1][0])/fft->vol;
  fft->four[1][1] = (fft->sclprim[2][2]*fft->sclprim[0][0] - fft->sclprim[0][2]*fft->sclprim[2][0])/fft->vol;
  fft->four[1][2] = (fft->sclprim[0][2]*fft->sclprim[1][0] - fft->sclprim[1][2]*fft->sclprim[0][0])/fft->vol;

  fft->four[2][0] = (fft->sclprim[1][0]*fft->sclprim[2][1] - fft->sclprim[2][0]*fft->sclprim[1][1])/fft->vol;
  fft->four[2][1] = (fft->sclprim[2][0]*fft->sclprim[0][1] - fft->sclprim[0][0]*fft->sclprim[2][1])/fft->vol;
  fft->four[2][2] = (fft->sclprim[0][0]*fft->sclprim[1][1] - fft->sclprim[1][0]*fft->sclprim[0][1])/fft->vol;

  fft->normfour[0] = sqrt(fft->four[0][0]*fft->four[0][0]+
		                      fft->four[0][1]*fft->four[0][1]+
		                      fft->four[0][2]*fft->four[0][2]);
  fft->normfour[1] = sqrt(fft->four[1][0]*fft->four[1][0]+
		                      fft->four[1][1]*fft->four[1][1]+
		                      fft->four[1][2]*fft->four[1][2]);
  fft->normfour[2] = sqrt(fft->four[2][0]*fft->four[2][0]+
		                      fft->four[2][1]*fft->four[2][1]+
		                      fft->four[2][2]*fft->four[2][2]);


  for (int i=0; i<dimension; i++) {
     for (int j=0; j<dimension; j++) {
       if(me == 0){
         if (logfile) fprintf(logfile,"vol = %lf, prim[%d][%d] = %lf, invprim[%d][%d] = %lf, sclprim[%d][%d] = %lf, four[%d][%d] = %lf\n",
                             fft->vol, i, j, fft->prim[i][j], i, j, fft->invprim[i][j], i, j, fft->sclprim[i][j], i, j, fft->four[i][j]);
         if (screen) fprintf(screen,"vol = %lf, prim[%d][%d] = %lf, invprim[%d][%d] = %lf, sclprim[%d][%d] = %lf, four[%d][%d] = %lf\n",
                             fft->vol, i, j, fft->prim[i][j], i, j, fft->invprim[i][j], i, j, fft->sclprim[i][j], i, j, fft->four[i][j]);
       }
     }
  }
}

void App::set_nonSchmid(int narg, char **arg)
{
  non_Schmid = 1;

  if (narg != 4) error->all(FLERR,"Illegal non_Schmid command");
  angle_to_110 = atof(arg[0]);
  w1 = atof(arg[1]);
  w2 = atof(arg[2]);
  w3 = atof(arg[3]);

  if(me == 0){
    if (logfile) fprintf(logfile,"angle_to_110= %lf, w1 = %lf, w2 = %lf, w3 = %lf\n", angle_to_110,w1,w2,w3);
    if (screen) fprintf(screen,"angle_to_110= %lf, w21 = %lf, w2 = %lf, w3 = %lf\n", angle_to_110,w1,w2,w3);
  }
}

/* ---------------------------------------------------------------------- */

void App::init()
{
  // error checks

  // if (solve == NULL && sweepflag == NOSWEEP)
  //   error->all(FLERR,"App needs a KMC or rejection KMC solver");
  // if (solve && sweepflag != NOSWEEP)
  //   error->all(FLERR,"App cannot use both a KMC and rejection KMC solver");

  // if (solve && allow_kmc == 0)
  //   error->all(FLERR,"KMC events are not implemented in app");
  // if (sweepflag != NOSWEEP && allow_rejection == 0)
  //   error->all(FLERR,"Rejection events are not implemented in app");
  // if (sweepflag != NOSWEEP && Lmask && allow_masking == 0)
  //   error->all(FLERR,"Mask logic not implemented in app");

  // if (nprocs > 1 && sectorflag == 0 && solve)
  //   error->all(FLERR,"Cannot use KMC solver in parallel with no sectors");
  // if (nprocs > 1 && sectorflag == 0 && sweepflag == RANDOM)
  //   error->all(FLERR,"Cannot use random rejection KMC in parallel with no sectors");
  // if (nprocs > 1 && sectorflag == 0 && sweepflag == RASTER)
  //   error->all(FLERR,"Cannot use raster rejection KMC in parallel with no sectors");
  // if (sectorflag && sweepflag == COLOR_STRICT)
  //   error->all(FLERR,"Cannot use color/strict rejection KMC with sectors");

  // if (sweepflag && dt_sweep == 0.0) error->all(FLERR,"App did not set dt_sweep");

  // setup RN generators, only on first init
  // ranapp is used for all options except sweep color/strict
  // setup ranapp so different on every proc
  // if color/strict, initialize per-lattice site seeds

  if (ranapp == NULL) {
    ranapp = new RandomPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    ranapp->reset(seed,me,100);
  }

  // fft-specific initialization, after general initialization
  // has to come before init_app()
  fft->init();
  // app-specific initialization, after general initialization
  init_app();

  // initialize comm, both for this proc's full fft and sectors
  // recall comm->init in case sectoring has changed

  if (comm == NULL) comm = new CommLattice(pfdd_p);
  comm->init(neighlayer,NULL);

  // initialize output

  output->init(time);
}

/* ---------------------------------------------------------------------- */

void App::setup()
{
  // app-specific setup, before propensities are computed

  setup_app();

  // initialize propensities for KMC solver within each set
  // comm insures ghost sites are up to date

  // if (solve) {
  //   comm->all();
  //   for (int i = 0; i < nset; i++) {
  //     for (int m = 0; m < set[i].nlocal; m++)
  // 	set[i].propensity[m] = site_propensity(set[i].site2i[m]);
  //     set[i].solve->init(set[i].nlocal,set[i].propensity);
  //   }
  // }

  // // convert per-sector time increment info to KMC params

  // if (sectorflag && solve) {
  //   if (tstop > 0.0) {
  //     Ladapt = false;
  //     dt_kmc = tstop;
  //   }

  //   if (nstop > 0.0) {
  //     Ladapt = true;
  //     double pmax = 0.0;
  //     for (int i = 0; i < nset; i++) {
  // 	int ntmp = set[i].solve->get_num_active();
  // 	if (ntmp > 0) {
  // 	  double ptmp = set[i].solve->get_total_propensity();
  // 	  ptmp /= ntmp;
  // 	  pmax = MAX(ptmp,pmax);
  // 	}
  //     }
  //     double pmaxall;
  //     MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
  //     if (pmaxall > 0.0) dt_kmc = nstop/pmaxall;
  //     else dt_kmc = stoptime-time;
  //   }

  //   dt_kmc = MIN(dt_kmc,stoptime-time);
  // }

  // // convert rejection info to rKMC params
  // // nloop and nselect are set whether sectoring is used or not
  // // if bothflag, list of active sets starts with nsector

  // if (sweepflag != NOSWEEP) {
  //   int first = 0;
  //   if (bothflag) first = nsector;

  //   if (sweepflag == RANDOM) {
  //     if (nstop > 0.0) {
  // 	for (int i = first; i < nset; i++) {
  // 	  set[i].nloop = 0;
  // 	  set[i].nselect = static_cast<int> (nstop*set[i].nlocal);
  // 	}
  //     }
  //     if (tstop > 0.0) {
  // 	double n = tstop / (dt_sweep/nglobal);
  // 	for (int i = first; i < nset; i++) {
  // 	  set[i].nloop = 0;
  // 	  set[i].nselect = static_cast<int> (n/nglobal * set[i].nlocal);
  // 	}
  //     }

  //   } else if (sweepflag == RASTER ||
  // 	       sweepflag == COLOR || sweepflag == COLOR_STRICT) {
  //     int n;
  //     if (nstop > 0.0) n = static_cast<int> (nstop);
  //     if (tstop > 0.0) n = static_cast<int> (tstop/dt_sweep);
  //     for (int i = first; i < nset; i++) {
  // 	set[i].nloop = n;
  // 	set[i].nselect = n * set[i].nlocal;
  //     }
  //   }

  //   double nme = 0.0;
  //   for (int i = first; i < nset; i++) nme += set[i].nselect;
  //   double ntotal;
  //   MPI_Allreduce(&nme,&ntotal,1,MPI_DOUBLE,MPI_SUM,world);

  //   dt_rkmc = ntotal/nglobal * dt_sweep;
  //   if (dt_rkmc == 0.0)
  //     error->all(FLERR,"Choice of sector stop led to no rKMC events");
  //   dt_rkmc = MIN(dt_rkmc,stoptime-time);
  // }

  // // setup sitelist if sweepflag = RANDOM
  // // do this every run since sector timestep could have changed

  // if (sweepflag == RANDOM) {
  //   memory->destroy(sitelist);
  //   int n = 0;
  //   for (int i = 0; i < nset; i++) n = MAX(n,set[i].nselect);
  //   memory->create(sitelist,n,"app:sitelist");
  // }

  // FFT setup
  fft->setup();

  // second stage of app-specific setup

  setup_end_app();

  // setup future output
  if (fft->me == 0) {
    if (screen) fprintf(screen,"Setting up output ...\n");
    if (logfile) fprintf(logfile,"Setting up output ...\n");
  }
  fft->temp_stor_data();
  fft->prep_forward();
  nextoutput = output->setup(time);
  fft->reset_data();
}

/* ---------------------------------------------------------------------- */

void App::iterate()
{
  timer->barrier_start(TIME_LOOP);

  solve->iterate();

  timer->barrier_stop(TIME_LOOP);
}

/* -----------------------------------------------------------------------
   Chooses external energy (schmid or non-Schmid) calculations
   ---------------------------------------------------------------------*/
void App::resolSS()
{
  if(non_Schmid == 0) fft->resolSS_Schmid();
  else fft->resolSS_non_Schmid();
}

/* -----------------------------------------------------------------------
   Chooses core energy
   ---------------------------------------------------------------------*/
void App::core_energy()
{
  if(core_flag == 1) fft->core_energy_perfect();
  else if(core_flag == 2) fft->core_energy_extended();
  else if(core_flag == 3) fft->core_energy_sine();
  else if(core_flag == 4) fft->core_energy_pyrII();// Claire changed 8/1/18: was core_energy_isf_usf
  else if(core_flag == 5) fft->core_energy_bcc_perfect();
  else if(core_flag == 6) fft->core_energy_USFE_angle_tau();  //added_HJ
  else if(core_flag == 7) fft->core_energy_mpea(); //added LTWS
}

/* ----------------------------------------------------------------------
   grow per-site arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void App::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  memory->grow(id,nmax,"app:id");
  memory->grow(xyz,nmax,3,"app:xyz");
  memory->grow(owner,nmax,"app:owner");
  memory->grow(index,nmax,"app:index");

  memory->grow(numneigh,nmax,"app:numneigh");
  if (maxneigh) memory->grow(neighbor,nmax,maxneigh,"app:neighbor");

  for (int i = 0; i < ninteger; i++)
    memory->grow(iarray[i],nmax,"app:iarray");
  for (int i = 0; i < ndouble; i++)
    memory->grow(darray[i],nmax,"app:darray");

  grow_app();
}

/* ----------------------------------------------------------------------
   add an owned site
   called from create_sites command
   grow arrays if necessary
 ------------------------------------------------------------------------- */

void App::add_site(tagint n, double x, double y, double z)
{
  if (nlocal == nmax) grow(0);

  id[nlocal] = n;
  xyz[nlocal][0] = x;
  xyz[nlocal][1] = y;
  xyz[nlocal][2] = z;

  owner[nlocal] = me;
  index[nlocal] = nlocal;

  for (int i = 0; i < ninteger; i++) iarray[i][nlocal] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][nlocal] = 0.0;

  nlocal++;
}

/* ----------------------------------------------------------------------
   add a ghost site
   called from create_sites command
   grow arrays if necessary
 ------------------------------------------------------------------------- */

void App::add_ghost(tagint n, double x, double y, double z,
			   int proc, int index_owner)
{
  if (nlocal+nghost == nmax) grow(0);

  int m = nlocal + nghost;

  id[m] = n;
  xyz[m][0] = x;
  xyz[m][1] = y;
  xyz[m][2] = z;

  owner[m] = proc;
  index[m] = index_owner;

  for (int i = 0; i < ninteger; i++) iarray[i][m] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][m] = 0.0;

  nghost++;
}

/* ----------------------------------------------------------------------
   set neighbor connectivity for owned site I
   nvalues = # of neighbors
   called from read_sites command
 ------------------------------------------------------------------------- */

void App::add_neighbors(int i, int nvalues, char **values)
{
  numneigh[i] = nvalues;
  for (int m = 0; m < nvalues; m++)
    neighbor[i][m] = atoi(values[m]);
}

/* ----------------------------------------------------------------------
   set values for owned site I
   called from read_sites command
 ------------------------------------------------------------------------- */

void App::add_values(int i, char **values)
{
  for (int m = 0; m < ninteger; m++) iarray[m][i] = atoi(values[m]);
  for (int m = 0; m < ndouble; m++) darray[m][i] = atof(values[m+ninteger]);
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   lo,hi = inclusive bounds
   5 possibilities:
     (1) i = i to i, (2) * = lo to hi,
     (3) i* = i to hi, (4) *j = lo to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void App::bounds(char *str, int lo, int hi, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = MAX(atoi(str),lo);
    nhi = MIN(atoi(str),hi);
  } else if (strlen(str) == 1) {
    nlo = lo;
    nhi = hi;
  } else if (ptr == str) {
    nlo = lo;
    nhi = MIN(atoi(ptr+1),hi);
  } else if (strlen(ptr+1) == 0) {
    nlo = MAX(atoi(str),lo);
    nhi = hi;
  } else {
    nlo = MAX(atoi(str),lo);
    nhi = MIN(atoi(ptr+1),hi);
  }
}

/* define initial dislocation structure */
void App::set_initial_sxtal(char *str)
{
  if (strcmp(str,"loop") == 0)
    setup_flag = 0;
  else if (strcmp(str,"none") == 0)
      setup_flag = 1;
  else if (strcmp(str,"edge") == 0)
      setup_flag = 2;
  else if (strcmp(str,"screw") == 0)
      setup_flag = 3;
  else if (strcmp(str,"mode1") == 0)
      setup_flag = 4;
  else if (strcmp(str,"1LOrthoFCC") == 0)
      setup_flag = 5;
  else if (strcmp(str,"1LNonOrthoFCC") == 0)
      setup_flag = 6;
  else if (strcmp(str,"1LInclinedFCC") == 0)
      setup_flag = 7;
  else if (strcmp(str,"1LOrthoBCC") == 0)
      setup_flag = 8;
  else if (strcmp(str,"2LOrtho") == 0)
      setup_flag = 9;
  else if (strcmp(str,"StraightOrtho") == 0)
      setup_flag = 10;
  else if (strcmp(str,"StraightNonOrtho") == 0)
      setup_flag = 11;
  else if (strcmp(str,"frank_read") == 0)
      setup_flag = 12;
  else if (strcmp(str,"mode2") == 0)
      setup_flag = 13;
  else if (strcmp(str,"3D3ShcpNotch") == 0)
      setup_flag = 14;
  else if (strcmp(str,"3D3ShcpBasal") == 0)
      setup_flag = 15;
  else
    error->all(FLERR,"Illegal setup command");
}

/* sets up dislocation configuration */
void App::initial_sxtal()
{
  if (setup_flag == 0)
    fft->initial_sxtal_loop();
  else if (setup_flag == 1)
    fft->initial_sxtal_NoLoop();
  else if (setup_flag == 2)
    fft->initial_sxtal_edge();
  else if (setup_flag == 3)
    fft->initial_sxtal_3SBCC();
  else if (setup_flag == 4)
    fft->initial_sxtal_mode1();
  else if (setup_flag == 5)
    fft->initial_sxtal_1LOrthoFCC();
  else if (setup_flag == 6)
    fft->initial_sxtal_1LNonOrthoFCC();
  else if (setup_flag == 7)
    fft->initial_sxtal_1LInclinedFCC();
  else if (setup_flag == 8)
    fft->initial_sxtal_1LOrthoBCC();
  else if (setup_flag == 9)
    fft->initial_sxtal_2LOrtho();
  else if (setup_flag == 10)
    fft->initial_sxtal_StraightOrtho();
  else if (setup_flag == 11)
    fft->initial_sxtal_StraightNonOrtho();
  else if (setup_flag == 12)
    fft->initial_sxtal_frank_read();
  else if (setup_flag == 13)
    fft->initial_sxtal_mode2();
  else if (setup_flag == 14)
    fft->initial_sxtal_3D3ShcpNotch();
  else if (setup_flag == 15)
    fft->initial_sxtal_3D3ShcpBasal();
  else
    error->all(FLERR,"Illegal setup command");
}

/* ----------------------------------------------------------------------
   print connectivity stats
 ------------------------------------------------------------------------- */

void App::print_connectivity()
{
  int i;

  tagint min = maxneigh;
  tagint max = 0;

  for (i = 0; i < nlocal; i++) {
    min = MIN(min,numneigh[i]);
    max = MAX(max,numneigh[i]);
  }

  tagint minall,maxall;
  MPI_Allreduce(&min,&minall,1,MPI_PFDD_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&max,&maxall,1,MPI_PFDD_TAGINT,MPI_MAX,world);

  tagint *count = new tagint[maxall+1];
  tagint *countall = new tagint[maxall+1];

  for (i = 0; i <= maxall; i++) count[i] = 0;

  for (i = 0; i < nlocal; i++) count[numneigh[i]]++;
  MPI_Allreduce(count,countall,maxall+1,MPI_PFDD_TAGINT,MPI_SUM,world);

  if (me == 0)
    for (i = minall; i <= maxall; i++) {
      if (screen)
	fprintf(screen,"  " TAGINT_FORMAT " sites have %d neighbors\n",\
		countall[i],i);
      if (logfile)
	fprintf(logfile,"  " TAGINT_FORMAT " sites have %d neighbors\n",
		countall[i],i);
    }

  delete [] count;
  delete [] countall;
}
