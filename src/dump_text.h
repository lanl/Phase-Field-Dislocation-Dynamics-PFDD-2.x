/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(text,DumpText)

#else

#ifndef PFDD_DUMP_TEXT_H
#define PFDD_DUMP_TEXT_H

#include "dump.h"

namespace PFDD_NS {

class DumpText : public Dump {
 public:
  DumpText(class PFDD_C *, int, char **);
  virtual ~DumpText();

 protected:
  int ioptional;             // where optional trailing args start

  int *vtype;                // type of each vector (INT, DOUBLE)
  int *vindex;               // index into int,double packs
  char **vformat;            // format string for each vector element

  int iregion;               // -1 if no region, else which region
  char *idregion;            // region ID

  int nthresh;               // # of defined threshholds
  int *thresh_array;         // array to threshold on for each nthresh
  int *thresh_op;            // threshold operation for each nthresh
  double *thresh_value;      // threshold value for each nthresh
  int *thresh_index;         // N index for iN and dN thresholds

  char *columns;             // text describing columns of dump output

  int nchoose;               // # of selected atoms
  int maxlocal;              // size of atom selection and variable arrays
  int *choose;               // lists of sites chosen for output
  double *dchoose;           // value for each atom to threshhold against
  int *clist;                // compressed list of indices of selected atoms

  int stress_flag;           // 1 if stress to be computed
  int strain_flag;           // 1 if strain to be computed
  int energy_flag;           // 1 if energy to be computed

  // private methods

  virtual void init_style();
  int count();
  void pack();
  void write_header(int, double);
  void write_data(int, double *);
  int parse_fields(int, char **);
  virtual int modify_param(int, char **);
  void compute();

  typedef void (DumpText::*FnPtrHeader)(int, double);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(int, double);
  void header_text(int, double);

  typedef void (DumpText::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_text(int, double *);

  typedef void (DumpText::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_xa(int);
  void pack_ya(int);
  void pack_za(int);
  void pack_i(int);
  void pack_j(int);
  void pack_k(int);
  void pack_theta(int); //Grad
  void pack_usfe1(int);
  void pack_usfe2(int);
  void pack_usfe3(int);
  void pack_delta(int); // projected order paramter
  void pack_ddelta(int); // gradient of delta
  void pack_nhn1(int);
  void pack_nhn2(int);
  void pack_picked_1(int);
  void pack_picked_2(int);
  void pack_dtkmc1(int);
  void pack_dtkmc2(int);
  void pack_gradx(int);
  void pack_grady(int);
  void pack_gradz(int);
  void pack_pxx(int);
  void pack_pyy(int);
  void pack_pzz(int);
  void pack_pxy(int);
  void pack_pxz(int);
  void pack_pyz(int);
  void pack_pxx_NoExt(int);
  void pack_pyy_NoExt(int);
  void pack_pzz_NoExt(int);
  void pack_pxy_NoExt(int);
  void pack_pxz_NoExt(int);
  void pack_pyz_NoExt(int);
  void pack_exx(int);
  void pack_eyy(int);
  void pack_ezz(int);
  void pack_exy(int);
  void pack_exz(int);
  void pack_eyz(int);

  void pack_xir1(int);
  void pack_xii1(int);
  void pack_xir2(int);
  void pack_xii2(int);
  void pack_xir3(int);
  void pack_xii3(int);
  void pack_xir4(int);
  void pack_xii4(int);
  void pack_xir5(int);
  void pack_xii5(int);
  void pack_xir6(int);
  void pack_xii6(int);
  void pack_xir7(int);
  void pack_xii7(int);
  void pack_xir8(int);
  void pack_xii8(int);
  void pack_xir9(int);
  void pack_xii9(int);
  void pack_xir10(int);
  void pack_xii10(int);
  void pack_xir11(int);
  void pack_xii11(int);
  void pack_xir12(int);
  void pack_xii12(int);
  void pack_xo1(int);
  void pack_xo2(int);
  void pack_xo3(int);
  void pack_xo4(int);
  void pack_xo5(int);
  void pack_xo6(int);
  void pack_dreg1(int);
  void pack_dreg2(int);
  void pack_dregp1_1(int);
  void pack_dregp1_2(int);
  void pack_dregp2_1(int);
  void pack_dregp2_2(int);
  void pack_int1(int);
  void pack_int2(int);
  void pack_int3(int);
  void pack_int4(int);
  void pack_int5(int);
  void pack_int6(int);
  void pack_ext1(int);
  void pack_ext2(int);
  void pack_ext3(int);
  void pack_ext4(int);
  void pack_ext5(int);
  void pack_ext6(int);
  void pack_lat1(int);
  void pack_lat2(int);
  void pack_lat3(int);
  void pack_lat4(int);
  void pack_lat5(int);
  void pack_lat6(int);
  void pack_data_real_1(int);
  void pack_ff1_xx(int);
  void pack_ff1_xy(int);
  void pack_ff1_xz(int);
  void pack_ff1_yy(int);
  void pack_ff1_yz(int);
  void pack_ff1_zz(int);
    void pack_usfe(int);

  void pack_energy(int);
  void pack_iarray(int);
  void pack_darray(int);

  void pack_iexx(int);
  void pack_ieyy(int);
  void pack_iezz(int);
  void pack_iexy(int);
  void pack_iexz(int);
  void pack_ieyz(int);

  void pack_apxx(int);
  void pack_apyy(int);
  void pack_apzz(int);
  void pack_apxy(int);
  void pack_apxz(int);
  void pack_apyz(int);
  void pack_tpxx(int);
  void pack_tpyy(int);
  void pack_tpzz(int);
  void pack_tpxy(int);
  void pack_tpxz(int);
  void pack_tpyz(int);
  void pack_ppxx(int);
  void pack_ppyy(int);
  void pack_ppzz(int);
  void pack_ppxy(int);
  void pack_ppxz(int);
  void pack_ppyz(int);
  void pack_upxx(int);
  void pack_upyy(int);
  void pack_upzz(int);
  void pack_upxy(int);
  void pack_upxz(int);
  void pack_upyz(int);
  void pack_rpxx(int);
  void pack_rpyy(int);
  void pack_rpzz(int);
  void pack_rpxy(int);
  void pack_rpxz(int);
  void pack_rpyz(int);
  void pack_phase(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid attribute in dump text command

UNDOCUMENTED

E: Dump requires propensity but no KMC solve performed

Only KMC solvers compute propensity for sites.

E: Region ID for dump text does not exist

UNDOCUMENTED

E: Dumping a quantity application does not support

The application defines what variables it supports.  You cannot
output a variable in a dump that isn't supported.

E: Invalid keyword in dump command

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Dump_modify region ID does not exist

UNDOCUMENTED

E: Threshold for a quantity application does not support

The application defines what variables it supports.  You cannot do a
threshold test with the dump command on a variable that isn't
supported.

E: Invalid dump_modify threshold operator

Self-explanatory.

*/
