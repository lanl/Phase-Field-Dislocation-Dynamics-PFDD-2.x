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
#include "Types.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "dump_text.h"
#include "app.h"
#include "fft.h"
#include "region.h"
#include "memory.h"
#include "error.h"
#include "output.h"

using namespace PFDD_NS;

// customize by adding keyword to 1st enum

enum{X,Y,Z,THETA,DELTAO,DDELTA,XIR1,XII1,XIR2,XII2,XIR3,XII3,XIR4,XI41,XIR5,XII5,XIR6,XII6,XIR7,XII7,XIR8,XII8,XIR9,XII9,XIR10,XII10,XIR11,XII11,XIR12,XII12,PXX,PYY,PZZ,PXY,PXZ,PYZ,EXX,EYY,EZZ,EXY,EXZ,EYZ,IARRAY,DARRAY};  // also in dump_image
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,TAGINT};           // also in dump_image

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

DumpText::DumpText(PFDD_C *pfdd_p, int narg, char **arg) : Dump(pfdd_p, narg, arg)
{
  // allow for default "x y z <xi>" output if no args specified

  int def = 0;
  char **argcopy;

  stress_flag = 1;
  strain_flag = 1;

  if (narg == 4) {
    def = 1;
    narg = 8;
    argcopy[4] = "x";
    argcopy[5] = "y";
    argcopy[6] = "z";
    argcopy[7] = new char[6];
    strcpy(argcopy[7],"<xi>");
  } else argcopy = arg;

  // size_one may be shrunk below if additional optional args exist

  size_one = narg - 4;
  vtype = new int[size_one];
  vindex = new int[size_one];
  vformat = new char*[size_one];
  pack_choice = new FnPtrPack[size_one];

  // process attributes
  // ioptional = start of additional optional args
  // only dump image style processes optional args

  ioptional = parse_fields(narg,argcopy);

  if (ioptional < narg && strcmp(style,"image") != 0)
    error->all(FLERR,"Invalid attribute in dump text command");
  size_one = ioptional - 4;

  // setup vformat strings, one per field

  for (int i = 0; i < size_one; i++) {
    char *format;
    if (vtype[i] == INT) format = "%d ";
    else if (vtype[i] == DOUBLE) format = "%g ";
    else if (vtype[i] == TAGINT) format = TAGINT_FORMAT " ";
    int n = strlen(format) + 1;
    vformat[i] = new char[n];
    strcpy(vformat[i],format);
  }

  // setup column string
  // change "site" to "type" to be LAMMPS compatible

  int n = 0;
  for (int iarg = 4; iarg < narg; iarg++) n += strlen(argcopy[iarg]) + 2;
  columns = new char[n];
  columns[0] = '\0';
  for (int iarg = 4; iarg < narg; iarg++) {
    if (strstr(argcopy[iarg],"site")) strcat(columns,"type");
    else strcat(columns,argcopy[iarg]);
    strcat(columns," ");
  }

  // delete argcopy if default output created

  if (def) {
    delete [] argcopy[7];
    delete [] argcopy;
  }

  // dump params

  iregion = -1;
  idregion = NULL;

  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;
  thresh_index = NULL;

  maxlocal = 0;
  choose = NULL;
  dchoose = NULL;
  clist = NULL;

  // setup function ptrs

  if (binary) header_choice = &DumpText::header_binary;
  else header_choice = &DumpText::header_text;

  if (binary) write_choice = &DumpText::write_binary;
  else write_choice = &DumpText::write_text;
}

/* ---------------------------------------------------------------------- */

DumpText::~DumpText()
{
  delete [] columns;

  delete [] idregion;
  memory->sfree(thresh_array);
  memory->sfree(thresh_op);
  memory->sfree(thresh_value);
  memory->sfree(thresh_index);

  memory->sfree(choose);
  memory->sfree(dchoose);
  memory->sfree(clist);

  delete [] vtype;
  delete [] vindex;
  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;
  delete [] pack_choice;
}

/* ---------------------------------------------------------------------- */

void DumpText::init_style()
{
  // error if propensity is dumped or used as threshold but doesn't exist
  // can't check until now, b/c input script may not have defined solver

  int flag = 0;
  // for (int i = 0; i < size_one; i++)
  //   if (pack_choice[i] == &DumpText::pack_propensity) flag = 1;
  // for (int i = 0; i < nthresh; i++)
  //   if (thresh_array[i] == PROPENSITY) flag = 1;
  // if (flag && !fft)
  //   error->all(FLERR,"Dump requires propensity but no KMC solve performed");

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = fft->find_region(idregion);
    if (iregion == -1) error->all(FLERR,"Region ID for dump text does not exist");
  }

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

int DumpText::count()
{
  int i;

  // grow choose arrays if needed

  int nlocal = app->nlocal;
  if (nlocal > maxlocal) {
    maxlocal = nlocal;

    memory->sfree(choose);
    memory->sfree(dchoose);
    memory->sfree(clist);
    choose = (int *) memory->smalloc(maxlocal*sizeof(int),"dump:choose");
    dchoose = (double *)
      memory->smalloc(maxlocal*sizeof(double),"dump:dchoose");
    clist = (int *) memory->smalloc(maxlocal*sizeof(int),"dump:clist");
  }

  // choose all local sites for output

  for (i = 0; i < nlocal; i++) choose[i] = 1;

  // un-choose if not in region

  if (iregion >= 0) {
    Region *region = fft->regions[iregion];
    double **xyz = app->xyz;
    for (i = 0; i < nlocal; i++)
      if (choose[i] && region->match(xyz[i][0],xyz[i][1],xyz[i][2]) == 0)
	choose[i] = 0;
  }

  // un-choose if any threshhold criterion isn't met

  if (nthresh) {
    double *ptr;
    double value;
    int nstride;

    for (int ithresh = 0; ithresh < nthresh; ithresh++) {

      if (thresh_array[ithresh] == X) {
	ptr = &app->xyz[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == Y) {
	ptr = &app->xyz[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == Z) {
	ptr = &app->xyz[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == IARRAY) {
	int index = thresh_index[ithresh];
	for (i = 0; i < nlocal; i++)
	  dchoose[i] = app->iarray[index][i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == DARRAY) {
	ptr = app->darray[thresh_index[ithresh]];
	nstride = 1;
      }

      // unselect sites that don't meet threshold criterion

      value = thresh_value[ithresh];

      if (thresh_op[ithresh] == LT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr >= value) choose[i] = 0;
      } else if (thresh_op[ithresh] == LE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr > value) choose[i] = 0;
      } else if (thresh_op[ithresh] == GT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr <= value) choose[i] = 0;
      } else if (thresh_op[ithresh] == GE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr < value) choose[i] = 0;
      } else if (thresh_op[ithresh] == EQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr != value) choose[i] = 0;
      } else if (thresh_op[ithresh] == NEQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr == value) choose[i] = 0;
      }
    }
  }

  // compress choose flags into clist
  // nchoose = # of selected atoms
  // clist[i] = local index of each selected atom

  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;

  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpText::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpText::compute()
{
  int ND = fft->dimension;
  for(int i=0;i<ND;i++){
    for(int j=0;j<ND;j++){
      fft->local_sigma[i][j] = 0.0;
      fft->avepst[i][j] = 0.0;
      fft->avepsts[i][j] = 0.0;
      fft->avepsd[i][j] = 0.0;
      fft->aveps[i][j] = 0.0;
      fft->ave_eps[i][j] = 0.0;
      fft->ave_sigma[i][j] = 0.0;
    }
  }
  if(stress_flag){
    fft->strain();
    fft->total_average_strain();
    fft->stress();
    output->comp_stress = 1;
    output->comp_strain = 1;
  }
  else if(strain_flag){
    fft->strain();
    fft->total_average_strain();
    output->comp_strain = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::write_header(int ndump, double time)
{
  if (multiproc) (this->*header_choice)(ndump,time);
  else if (me == 0) (this->*header_choice)(ndump,time);
}

/* ---------------------------------------------------------------------- */

void DumpText::write_data(int n, double *buf)
{
  (this->*write_choice)(n,buf);
}

/* ---------------------------------------------------------------------- */

void DumpText::header_binary(int ndump, double time)
{
  fwrite(&idump,sizeof(int),1,fp);
  fwrite(&time,sizeof(double),1,fp);
  fwrite(&ndump,sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) {
    int one = 1;
    fwrite(&one,sizeof(int),1,fp);
  } else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpText::header_text(int ndump, double time)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d %10g\n",idump,time);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",app->sa*boxxlo,app->sa*boxxhi);
  fprintf(fp,"%g %g\n",app->sb*boxylo,app->sb*boxyhi);
  fprintf(fp,"%g %g\n",app->sc*boxzlo,app->sc*boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpText::write_binary(int n, double *buf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(buf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpText::write_text(int n, double *buf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT)
	fprintf(fp,vformat[j],static_cast<int> (buf[m]));
      else if (vtype[j] == DOUBLE)
	fprintf(fp,vformat[j],buf[m]);
      else if (vtype[j] == TAGINT)
	fprintf(fp,vformat[j],static_cast<tagint> (buf[m]));
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpText::parse_fields(int narg, char **arg)
{
  int i;
  for (int iarg = 4; iarg < narg; iarg++) {
    i = iarg-4;

    if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &DumpText::pack_x;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &DumpText::pack_y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &DumpText::pack_z;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"theta") == 0) {  //Grad
      pack_choice[i] = &DumpText::pack_theta;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"usfe1") == 0) {
      pack_choice[i] = &DumpText::pack_usfe1;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"usfe2") == 0) {
      pack_choice[i] = &DumpText::pack_usfe2;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"usfe3") == 0) {
      pack_choice[i] = &DumpText::pack_usfe3;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"deltao") == 0) {  //projected
      pack_choice[i] = &DumpText::pack_delta;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"ddelta") == 0) {  //gradient of ddelta
      pack_choice[i] = &DumpText::pack_ddelta;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"pxx") == 0) {
      pack_choice[i] = &DumpText::pack_pxx;
      vtype[i] = DOUBLE;
      stress_flag = 1;
    } else if (strcmp(arg[iarg],"pyy") == 0) {
      pack_choice[i] = &DumpText::pack_pyy;
      vtype[i] = DOUBLE;
      stress_flag = 1;
    } else if (strcmp(arg[iarg],"pzz") == 0) {
      pack_choice[i] = &DumpText::pack_pzz;
      vtype[i] = DOUBLE;
      stress_flag = 1;
    } else if (strcmp(arg[iarg],"pxy") == 0) {
      pack_choice[i] = &DumpText::pack_pxy;
      vtype[i] = DOUBLE;
      stress_flag = 1;
    } else if (strcmp(arg[iarg],"pxz") == 0) {
      pack_choice[i] = &DumpText::pack_pxz;
      vtype[i] = DOUBLE;
      stress_flag = 1;
    } else if (strcmp(arg[iarg],"pyz") == 0) {
      pack_choice[i] = &DumpText::pack_pyz;
      vtype[i] = DOUBLE;
      stress_flag = 1;
    } else if (strcmp(arg[iarg],"exx") == 0) {
      pack_choice[i] = &DumpText::pack_exx;
      vtype[i] = DOUBLE;
      strain_flag = 1;
    } else if (strcmp(arg[iarg],"eyy") == 0) {
      pack_choice[i] = &DumpText::pack_eyy;
      vtype[i] = DOUBLE;
      strain_flag = 1;
    } else if (strcmp(arg[iarg],"ezz") == 0) {
      pack_choice[i] = &DumpText::pack_ezz;
      vtype[i] = DOUBLE;
      strain_flag = 1;
    } else if (strcmp(arg[iarg],"exy") == 0) {
      pack_choice[i] = &DumpText::pack_exy;
      vtype[i] = DOUBLE;
      strain_flag = 1;
    } else if (strcmp(arg[iarg],"exz") == 0) {
      pack_choice[i] = &DumpText::pack_exz;
      vtype[i] = DOUBLE;
      strain_flag = 1;
    } else if (strcmp(arg[iarg],"eyz") == 0) {
      pack_choice[i] = &DumpText::pack_eyz;
      vtype[i] = DOUBLE;
      strain_flag = 1;
    } else if (strcmp(arg[iarg],"xir1") == 0) {
      pack_choice[i] = &DumpText::pack_xir1;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii1") == 0) {
      pack_choice[i] = &DumpText::pack_xii1;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir2") == 0) {
      pack_choice[i] = &DumpText::pack_xir2;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii2") == 0) {
      pack_choice[i] = &DumpText::pack_xii2;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir3") == 0) {
      pack_choice[i] = &DumpText::pack_xir3;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii3") == 0) {
      pack_choice[i] = &DumpText::pack_xii3;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir4") == 0) {
      pack_choice[i] = &DumpText::pack_xir4;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii4") == 0) {
      pack_choice[i] = &DumpText::pack_xii4;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir5") == 0) {
      pack_choice[i] = &DumpText::pack_xir5;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii5") == 0) {
      pack_choice[i] = &DumpText::pack_xii5;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir6") == 0) {
      pack_choice[i] = &DumpText::pack_xir6;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii6") == 0) {
      pack_choice[i] = &DumpText::pack_xii6;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir7") == 0) {
      pack_choice[i] = &DumpText::pack_xir7;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii7") == 0) {
      pack_choice[i] = &DumpText::pack_xii7;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir8") == 0) {
      pack_choice[i] = &DumpText::pack_xir8;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii8") == 0) {
      pack_choice[i] = &DumpText::pack_xii8;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir9") == 0) {
      pack_choice[i] = &DumpText::pack_xir9;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii9") == 0) {
      pack_choice[i] = &DumpText::pack_xii9;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir10") == 0) {
      pack_choice[i] = &DumpText::pack_xir10;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii10") == 0) {
      pack_choice[i] = &DumpText::pack_xii10;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir11") == 0) {
      pack_choice[i] = &DumpText::pack_xir11;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii11") == 0) {
      pack_choice[i] = &DumpText::pack_xii11;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xir12") == 0) {
      pack_choice[i] = &DumpText::pack_xir12;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xii12") == 0) {
      pack_choice[i] = &DumpText::pack_xii12;
      vtype[i] = DOUBLE;

    // integer value = iN
    // double value = dN

    } else if (arg[iarg][0] == 'i') {
      pack_choice[i] = &DumpText::pack_iarray;
      vtype[i] = INT;
      vindex[i] = atoi(&arg[iarg][1]);
      if (latticeflag && (vindex[i] < 1 || vindex[i] > app->ninteger))
	error->all(FLERR,"Invalid keyword in dump command");
      vindex[i]--;
    } else if (arg[iarg][0] == 'd') {
      pack_choice[i] = &DumpText::pack_darray;
      vtype[i] = DOUBLE;
      vindex[i] = atoi(&arg[iarg][1]);
      if (latticeflag && (vindex[i] < 1 || vindex[i] > app->ndouble))
	error->all(FLERR,"Invalid keyword in dump command");
      vindex[i]--;

    } else return iarg;
  }

  return narg;
}

/* ---------------------------------------------------------------------- */

int DumpText::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) iregion = -1;
    else {
      iregion = fft->find_region(arg[1]);
      if (iregion == -1) error->all(FLERR,"Dump_modify region ID does not exist");
      int n = strlen(arg[1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[1]);
    }
    return 2;

  } else if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
	memory->sfree(thresh_array);
	memory->sfree(thresh_op);
	memory->sfree(thresh_value);
	memory->sfree(thresh_index);
	thresh_array = NULL;
	thresh_op = NULL;
	thresh_value = NULL;
	thresh_index = NULL;
      }
      nthresh = 0;
      return 2;
    }

    if (narg < 4) error->all(FLERR,"Illegal dump_modify command");

    // grow threshold arrays

    thresh_array = (int *)
      memory->srealloc(thresh_array,(nthresh+1)*sizeof(int),
		       "dump:thresh_array");
    thresh_op = (int *)
      memory->srealloc(thresh_op,(nthresh+1)*sizeof(int),
		       "dump:thresh_op");
    thresh_value = (double *)
      memory->srealloc(thresh_value,(nthresh+1)*sizeof(double),
		       "dump:thresh_value");
    thresh_index = (int *)
      memory->srealloc(thresh_index,(nthresh+1)*sizeof(int),
		       "dump:thresh_index");

    // set keyword type of threshold
    // customize by adding to if statement

    if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
    else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
    else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;

    // TO DO
    // Add order parameter, stress and strain


    // integer value = iN
    // double value = dN

    else if (arg[1][0] == 'i') {
      thresh_array[nthresh] = IARRAY;
      thresh_index[nthresh] = atoi(&arg[1][1]);
      if (thresh_index[nthresh] < 1 ||
	  thresh_index[nthresh] > app->ninteger)
	error->all(FLERR,"Threshold for a quantity application does not support");
      thresh_index[nthresh]--;
    } else if (arg[1][0] == 'd') {
      thresh_array[nthresh] = DARRAY;
      thresh_index[nthresh] = atoi(&arg[1][1]);
      if (thresh_index[nthresh] < 1 ||
	  thresh_index[nthresh] > app->ndouble)
	error->all(FLERR,"Threshold for a quantity application does not support");
      thresh_index[nthresh]--;

    } else error->all(FLERR,"Invalid dump_modify threshold operator");

    // set operation type of threshold

    if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
    else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
    else error->all(FLERR,"Invalid dump_modify threshold operator");

    // set threshold value

    thresh_value[nthresh] = atof(arg[3]);

    nthresh++;
    return 4;
  }

  return 0;
}


/* ----------------------------------------------------------------------
   one method for every keyword dump can output
   the site quantity is packed into buf starting at n with stride size_one
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */

void DumpText::pack_x(int n)
{
  double **xyz = app->xyz;
  // double norm = sqrt(fft->four[0][0]*fft->four[0][0]+
  // 		     fft->four[0][1]*fft->four[0][1]+
  // 		     fft->four[0][2]*fft->four[0][2]);

  for (int i = 0; i < nchoose; i++) {
    // buf[n] = (fft->four[0][0]*xyz[clist[i]][0]+
    // 	      fft->four[0][1]*xyz[clist[i]][1]+
    // 	      fft->four[0][2]*xyz[clist[i]][2])/norm;
    buf[n] = (fft->sclprim[0][0]*xyz[clist[i]][0]+
	      fft->sclprim[0][1]*xyz[clist[i]][1]+
	      fft->sclprim[0][2]*xyz[clist[i]][2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_y(int n)
{
  double **xyz = app->xyz;
  // double norm = sqrt(fft->four[1][0]*fft->four[1][0]+
  // 		     fft->four[1][1]*fft->four[1][1]+
  // 		     fft->four[1][2]*fft->four[1][2]);

  for (int i = 0; i < nchoose; i++) {
    // buf[n] = (fft->four[1][0]*xyz[clist[i]][0]+
    // 	      fft->four[1][1]*xyz[clist[i]][1]+
    // 	      fft->four[1][2]*xyz[clist[i]][2])/norm;
    buf[n] = (fft->sclprim[1][0]*xyz[clist[i]][0]+
	      fft->sclprim[1][1]*xyz[clist[i]][1]+
	      fft->sclprim[1][2]*xyz[clist[i]][2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_z(int n)
{
  double **xyz = app->xyz;
  // double norm = sqrt(fft->four[2][0]*fft->four[2][0]+
  // 		     fft->four[2][1]*fft->four[2][1]+
  // 		     fft->four[2][2]*fft->four[2][2]);

  for (int i = 0; i < nchoose; i++) {
    // buf[n] = (fft->four[2][0]*xyz[clist[i]][0]+
    // 	      fft->four[2][1]*xyz[clist[i]][1]+
    // 	      fft->four[2][2]*xyz[clist[i]][2])/norm;
    buf[n] = (fft->sclprim[2][0]*xyz[clist[i]][0]+
	      fft->sclprim[2][1]*xyz[clist[i]][1]+
	      fft->sclprim[2][2]*xyz[clist[i]][2]);
    n += size_one;
  }
}
/* ---------------------------------------------------------------------- */

void DumpText::pack_theta(int n)
{
  for (int i = 0; i < nchoose; i++) {
    //buf[n] = fft->theta[clist[i]];
    buf[n] = app->compute_theta(clist[i],1);
    n += size_one;
  }
}

void DumpText::pack_usfe1(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_usfe(clist[i],1);
    n += size_one;
  }
}

void DumpText::pack_usfe2(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_usfe(clist[i],2);
    n += size_one;
  }
}

void DumpText::pack_usfe3(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_usfe(clist[i],3);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_delta(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_delta(clist[i],1);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_ddelta(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_ddelta(clist[i],1);
    n += size_one;
  }
}


/* ---------------------------------------------------------------------- */

void DumpText::pack_xir1(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],1);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii1(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],2);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir2(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],3);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii2(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],4);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir3(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],5);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii3(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],6);
    n += size_one;
  }
}
/* ---------------------------------------------------------------------- */

void DumpText::pack_xir4(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],7);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii4(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],8);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir5(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],9);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii5(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],10);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir6(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],11);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii6(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],12);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir7(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],13);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii7(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],14);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir8(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],15);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii8(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],16);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir9(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],17);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii9(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],18);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir10(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],19);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii10(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],20);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir11(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],21);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii11(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],22);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xir12(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],23);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_xii12(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_mean_xi(clist[i],24);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_pxx(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_stress(clist[i],1);
    n += size_one;
  }
}
/* ---------------------------------------------------------------------- */
void DumpText::pack_pyy(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_stress(clist[i],2);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */
void DumpText::pack_pzz(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_stress(clist[i],3);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */
void DumpText::pack_pxy(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_stress(clist[i],4);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */
void DumpText::pack_pxz(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_stress(clist[i],5);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */
void DumpText::pack_pyz(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_stress(clist[i],6);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_exx(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_strain(clist[i],1);
    n += size_one;
  }
}
/* ---------------------------------------------------------------------- */
void DumpText::pack_eyy(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_strain(clist[i],2);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */
void DumpText::pack_ezz(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_strain(clist[i],3);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */
void DumpText::pack_exy(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_strain(clist[i],4);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */
void DumpText::pack_exz(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_strain(clist[i],5);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */
void DumpText::pack_eyz(int n)
{
  for (int i = 0; i < nchoose; i++) {
    buf[n] = app->compute_strain(clist[i],6);
    n += size_one;
  }
}
/* ---------------------------------------------------------------------- */

void DumpText::pack_iarray(int n)
{
  int *ivec = app->iarray[vindex[n]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = ivec[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_darray(int n)
{
  double *dvec = app->darray[vindex[n]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = dvec[clist[i]];
    n += size_one;
  }
}
