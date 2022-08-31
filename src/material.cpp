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

#include "Types.h"
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "material.h"
#include "fft.h"
#include "memory.h"
#include "error.h"

#define WORD 256

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

Material::Material(PFDD_C *pfdd_p, int narg, char **arg, int N1, int N2, int N3, int NS) : Pointers(pfdd_p)
{
  // parse style arg
  style = new char [WORD];
  strcpy(style,"user");

  if (narg < 2) error->all(FLERR,"Illegal material command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"a") == 0) {
      a = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"mu") == 0) {
      mu = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"young") == 0) {
      young = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"c0") == 0) {
      c0 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"c1") == 0) {
      c1 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"c2") == 0) {
      c2 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"c3") == 0) {
      c3 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"c4") == 0) {
      c4 = atof(arg[iarg+2]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"a1") == 0) {
      a1 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"a3") == 0) {
      a3 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"isf") == 0) {
      isf = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"usf") == 0) {
      usf = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"An") == 0) {
      An = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"Cn") == 0) {
      Cn = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"preset") == 0) {
      set_constants(arg[iarg+1]);
      iarg+=2;
      strcpy(style,"preset");
    }
    else if (strcmp(arg[iarg],"mpea") == 0){
      read_mpea_file(arg[iarg+1], N1, N2, N3, NS);
      iarg+=2;
    }
  }

  C44 = mu; //Pa
  nu = young/2.0/(mu)-1.0;
  C12 = 2.0*nu*C44/(1.0-2.0*nu); //Pa
  C11 = 2.0*C44+C12; //Pa
  ll = C12;

  S11 = 1.0/(young);
  S12 = -nu/(young);
  S44 = 2*(S11-S12);

  // moved b into the preset material cases so it work for both BCC and FCC
  // b = sqrt(3.0)*(a)/2.0; //meters-->a is in meters
  // b = (a)/sqrt(2.0);
  dslip = 0.8165; /*dslip=d/b and d=2(sqrt(3)a/2) or d = Constant*b*/
}

/* ---------------------------------------------------------------------- */

Material::~Material()
{
  //memory->destroy(basis);
}

/* ----------------------------------------------------------------------
   add a basis atom to list
   x,y,z = fractional coords within unit cell
------------------------------------------------------------------------- */

void Material::set_constants(const char *style)
{
 if (strcmp(style,"Ni") == 0) {
   a = 3.52E-10;
   b = (a)/sqrt(2.0);
   mu = 99.7E9;
   young = 255.7E9;

   c0 = 503.144E-3;
   c1 = -106.533E-3;
   c2 = -83.741E-3;
   c3 = 28.777E-3;
   c4 = -2.916E-3;
   a1 = -123.917E-3;
   a3 = -38.533E-3;

   isf = 144.5E-3;  //25.0E-3;
   usf = 289.0E-3; //270.0E-3;

 }
 else if (strcmp(style,"Al") == 0) {
   a = 4.04E-10;
   b = (a)/sqrt(2.0);
   mu = 28.4E9;
   young = 76.0E9;

   c0 = 149.56346E-3;
   c1 = -4.23464E-3;
   c2 = -74.39331E-3;
   c3 = 33.84682E-3;
   c4 = -2.39438E-3;
   a1 = 67.75268E-3;
   a3 = -30.88995E-3;

   isf = 140.2E-3;
   usf = 177.0E-3;
 }
 else if (strcmp(style,"Au") == 0) {

   a = 4.17E-10;
   b = (a)/sqrt(2.0);
   mu = 18.8E9;
   young = 54.0E9;

   c0 = 161.832E-3;
   c1 = -40.659E-3;
   c2 = -17.320E-3;
   c3 = 3.197E-3;
   c4 = 0.472E-3;
   a1 = -64.718E-3;
   a3 = -13.282E-3;

   isf = 27.9E-3; //37.0E-3;
   usf = 66.5E-3;
 }
 else if (strcmp(style,"Cu") == 0) {
   a = 3.64E-10; //meters
   b = (a)/sqrt(2.0);
   mu = 54.5E9; //Pascal
   young = 144.4E9; //Pascal

   c0 = 344.948E-3;
   c1 = -85.399E-3;
   c2 = -40.421E-3;
   c3 = 12.396E-3;
   c4 = -0.674E-3;
   a1 = -131.393E-3;
   a3 = -19.816E-3;

   isf = 38.5E-3;
   usf = 163.7E-3;
 }
 else if (strcmp(style,"Pd") == 0) {
   a = 3.950E-10;
   b = (a)/sqrt(2.0);
   mu = 45.7E9;
   young = 125.8E9;

   c0 = 308.65990E-3;
   c1 = -63.66363E-3;
   c2 = -48.76003E-3;
   c3 = 12.024463E-3;
   c4 = -1.11869E-3;
   a1 = -61.34161E-3;
   a3 = -23.33224E-3;

   isf = 138.1E-3;
   usf = 197.9E-3;
 }
 else if (strcmp(style,"Ag") == 0) {
   a = 4.16E-10;
   b = (a)/sqrt(2.0);
   mu = 27.5E9;
   young = 74.8E9;

   c0 = 179.573E-3;
   c1 = -32.837E-3;
   c2 = -39.578E-3;
   c3 = 14.700E-3;
   c4 = -0.989E-3;
   a1 = -43.791E-3;
   a3 = -17.254E-3;

   isf = 17.8E-3; //37.0E-3;
   usf = 100.4E-3;
 }
 else if (strcmp(style,"Pt") == 0) {
   a = 3.98E-10;
   b = (a)/sqrt(2.0);
   mu = 53.2E9;
   young = 149.0E9;

   c0 = 223.27233E-3;
   c1 = -15.38077E-3;
   c2 = -75.36747E-3;
   c3 = 21.70531E-3;
   c4 = -2.52876E-3;
   a1 = 58.60391E-3;
   a3 = -43.19019E-3;

   isf = 253.7E-3;
   usf = 257.5E-3;
 }
 else if (strcmp(style,"Rh") == 0) {
   a = 3.84E-10;
   b = (a)/sqrt(2.0);
   mu = 154.1E9;
   young = 384.9E9;

   c0 = 556.35590E-3;
   c1 = -46.42878E-3;
   c2 = -183.45940E-3;
   c3 = 55.25737E-3;
   c4 = -5.04421E-3;
   a1 = 1.44561E-3;
   a3 = -72.22364E-3;

   isf = 184.6E-3;
   usf = 454.1E-3;
 }
 else if (strcmp(style,"Ir") == 0) {
   a = 3.88E-10;
   b = (a)/sqrt(2.0);
   mu = 222.3E9;
   young = 549.7E9;

   c0 = 617.61574E-3;
   c1 = 3.28767E-3;
   c2 = -283.1158E-3;
   c3 = 96.12464E-3;
   c4 = -10.56298E-3;
   a1 = 157.20734E-3;
   a3 = -109.4680E-3;

   isf = 324.4E-3;
   usf = 614.9E-3;
 }
 else if (strcmp(style,"Mg") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 112.6205E-3;
   c1 = 0.7200E-3;
   c2 = -58.3781E-3;
   c3 = 28.1313E-3;
   c4 = -3.9212E-3;
   a1 = 25.2325E-3;
   a3 = -23.8545E-3;

   isf = 29.3216E-3;
   usf = 88.4906E-3;

   aa0 = 216.6E-3; //283.0E-3;
   aa1 = -82.84E-3; //-151.7E-3;
   aa2 = -116.8E-3; //-97.77E-3;
   aa3 = -5.873E-3; //-23.28E-3;
   aa4 = -11.27E-3; //-7.528E-3;
   bb1 = -51.94E-3; //-110.0E-3;
   bb2 = 22.04E-3; //43.44E-3;
   bb3 = 9.44E-3; //2.61E-3;
   bb4 = -0.2501E-3; //2.014E-3;
 }
 else if (strcmp(style,"MgY") == 0) { /*Magnesium47-Yttrium Basal */
     //Parameters from Pei paper: Ab initio and atomistic study of generalized stacking fault energies in Mg and Mg47-Y alloys.
   a = 3.21E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 126.3042E-3;
   c1 = -21.6452E-3;
   c2 = -29.3248E-3;
   c3 = 11.4053E-3;
   c4 = -1.2093E-3;
   a1 = -20.2756E-3;
   a3 = -10.9327E-3;

   isf = 30.0000E-3;
   usf = 83.0000E-3;
 }
 else if (strcmp(style,"Ti") == 0) {
   a = 2.95E-10;
   // need to add Burges vector
   // b = * a
   mu = 44.6E9;
   young = 118.67E9;

   c0 = 112.6205E-3; /* Basal constants for Mg, don't have values for Ti yet */
   c1 = 0.7200E-3;
   c2 = -58.3781E-3;
   c3 = 28.1313E-3;
   c4 = -3.9212E-3;
   a1 = 25.2325E-3;
   a3 = -23.8545E-3;

   isf = 333.565E-3;
   usf = 622.5063E-3;

   aa0 = 427.0E-3;
   aa1 = -165.9E-3;
   aa2 = -221.9E-3;
   aa3 = -26.29E-3;
   aa4 = -13.42E-3;
   bb1 = -53.09E-3;
   bb2 = 119.1E-3;
   bb3 = -14.25E-3;
   bb4 = -17.77E-3;
 }
/* The following MgX alloys ( where X = Al,Ca,Ce,Gd,Li,Si,Sn,Zn,Zr) are taken from A. Moitra, S. G. Kim, and M. F. Horstemeyer, “Solute effect on basal and prismatic slip systems of Mg,” J. Phys. Condens. Matter, vol. 26, no. 44, 2014. Constants a, mu, and young are the same as for Mg so we are only comparing the effects of GSFE differences on SFW */
 else if (strcmp(style,"MgAl") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 265.564638E-3;
   c1 = -96.111676E-3;
   c2 = 11.786691E-3;
   c3 = -7.169259E-3;
   c4 = 1.507026E-3;
   a1 = -161.871180E-3;
   a3 = 0.742991E-3;

   isf = 23.52941E-3;
   usf = 83.83E-3;

   aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
   aa1 = -82.84E-3;
   aa2 = -116.8E-3;
   aa3 = -5.873E-3;
   aa4 = -11.27E-3;
   bb1 = -51.94E-3;
   bb2 = 22.04E-3;
   bb3 = 9.44E-3;
   bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"MgCa") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 310.853100E-3;
   c1 = -114.443600E-3;
   c2 = 16.958600E-3;
   c3 = -11.221200E-3;
   c4 = 2.570400E-3;
   a1 = -191.432500E-3;
   a3 = 2.798000E-3;

   isf = 31.56863E-3;
   usf = 100.39216E-3;

   aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
   aa1 = -82.84E-3;
   aa2 = -116.8E-3;
   aa3 = -5.873E-3;
   aa4 = -11.27E-3;
   bb1 = -51.94E-3;
   bb2 = 22.04E-3;
   bb3 = 9.44E-3;
   bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"MgCe") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 241.210583E-3;
   c1 = -72.302347E-3;
   c2 = -13.746107E-3;
   c3 = 4.705126E-3;
   c4 = 0.523234E-3;
   a1 = 0.523234E-3;
   a3 = -7.398787E-3;

   isf = 1.17647E-3;
   usf = 95.19608E-3;

     aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
     aa1 = -82.84E-3;
     aa2 = -116.8E-3;
     aa3 = -5.873E-3;
     aa4 = -11.27E-3;
     bb1 = -51.94E-3;
     bb2 = 22.04E-3;
     bb3 = 9.44E-3;
     bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"MgGd") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 285.763113E-3;
   c1 = -95.960841E-3;
   c2 = 0.671163E-3;
   c3 = -1.588154E-3;
   c4 = 0.848439E-3;
   a1 = -159.569265E-3;
   a3 = -1.518973E-3;

   isf = 11.66667E-3;
   usf = 102.64706E-3;

     aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
     aa1 = -82.84E-3;
     aa2 = -116.8E-3;
     aa3 = -5.873E-3;
     aa4 = -11.27E-3;
     bb1 = -51.94E-3;
     bb2 = 22.04E-3;
     bb3 = 9.44E-3;
     bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"MgLi") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 281.067152E-3;
   c1 = -103.882445E-3;
   c2 = 16.615828E-3;
   c3 = -9.525039E-3;
   c4 = 1.569110E-3;
   a1 = -166.140969E-3;
   a3 = 4.570589E-3;

   isf = 49.01961E-3;
   usf = 99.41176E-3;

     aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
     aa1 = -82.84E-3;
     aa2 = -116.8E-3;
     aa3 = -5.873E-3;
     aa4 = -11.27E-3;
     bb1 = -51.94E-3;
     bb2 = 22.04E-3;
     bb3 = 9.44E-3;
     bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"MgSi") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 237.616625E-3;
   c1 = -85.750320E-3;
   c2 = 10.157547E-3;
   c3 = -5.707935E-3;
   c4 = 1.062990E-3;
   a1 = -148.577990E-3;
   a3 = -1.072479E-3;

   isf = 14.90196E-3;
   usf = 66.86275E-3;

     aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
     aa1 = -82.84E-3;
     aa2 = -116.8E-3;
     aa3 = -5.873E-3;
     aa4 = -11.27E-3;
     bb1 = -51.94E-3;
     bb2 = 22.04E-3;
     bb3 = 9.44E-3;
     bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"MgSn") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 265.383735E-3;
   c1 = -96.234971E-3;
   c2 = 12.703843E-3;
   c3 = -6.926203E-3;
   c4 = 1.011138E-3;
   a1 = -168.851904E-3;
   a3 = -0.639299E-3;

   isf = 12.25490E-3;
   usf = 72.84314E-3;

     aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
     aa1 = -82.84E-3;
     aa2 = -116.8E-3;
     aa3 = -5.873E-3;
     aa4 = -11.27E-3;
     bb1 = -51.94E-3;
     bb2 = 22.04E-3;
     bb3 = 9.44E-3;
     bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"MgZn") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 239.664844E-3;
   c1 = -83.074958E-3;
   c2 = 4.673015E-3;
   c3 = -3.589689E-3;
   c4 = 1.080065E-3;
   a1 = -132.943686E-3;
   a3 = -1.228621E-3;

   isf = 33.03922E-3;
   usf = 85.09804E-3;

     aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
     aa1 = -82.84E-3;
     aa2 = -116.8E-3;
     aa3 = -5.873E-3;
     aa4 = -11.27E-3;
     bb1 = -51.94E-3;
     bb2 = 22.04E-3;
     bb3 = 9.44E-3;
     bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"MgZr") == 0) {
   a = 3.19E-10;
   // need to add Burges vector
   // b = * a
   mu = 18.0E9;
   young = 46.8E9;

   c0 = 288.840762E-3;
   c1 = -104.774131E-3;
   c2 = 13.194615E-3;
   c3 = -8.018332E-3;
   c4 = 1.683542E-3;
   a1 = -168.501923E-3;
   a3 = 2.901528E-3;

   isf = 40.00000E-3;
   usf = 100.78431E-3;

     aa0 = 216.6E-3; /* constants for Mg Pyr2, need to be added for MgX */
     aa1 = -82.84E-3;
     aa2 = -116.8E-3;
     aa3 = -5.873E-3;
     aa4 = -11.27E-3;
     bb1 = -51.94E-3;
     bb2 = 22.04E-3;
     bb3 = 9.44E-3;
     bb4 = -0.2501E-3;
 }
 else if (strcmp(style,"Ta") == 0) {
   a = 3.3058E-10;
   b = sqrt(3.0)*(a)/2.0;
   mu = 70.67E9;
   young = 189.24E9;
   usf = 950.21E-3;

   //  Estimate the Ta USFE (MRSSP angle, tau)--SNAP
   //fitted angle: -30^o to 30^o, stresses 0.25, 0.5, 0.75GPa
   a_slope = 51.59816407 ;
   a_b = -2.27173784 ;
   b_b = 0.9979452 ;
   c_slope = -4.16536 ;
   c_b = 1133.961364 ;
 }

 else if (strcmp(style,"Fe") == 0) {
  a = 2.8665E-10;
  b = sqrt(3.0)*(a)/2.0;
  mu = 87.2796E9;
  young = 225.1814E9;
  usf = 744.14E-3;

  // Estimate the Fe USFE (MRSSP angle, tau)-- EAM, Proville
  //fitted angle: -30^o to 30^o, stresses 0.1, 0.25, 0.4 GPa
  // Estimate the Fe USFE (MRSSP angle, tau)
  a_slope = 15.524609942496577 ;
  a_b = 0.004388435835439675;
  b_b = 0.9934680885709689;
  c_slope = -0.9065590784654822;
  c_b = 944.0633238797747  ;
  }
  else if (strcmp(style,"Nb") == 0) {
    a = 3.324E-10;
    b = 2.859E-10;
    mu = 39.64E9;
    young = 110.33E9;
    usf = 676.775E-3;
    An = usf;
  }
  else if (strcmp(style,"Mo") ==  0){
    a=3.16E-10;
    b= sqrt(3.0)*a/2.0;
    mu = 119.29E9;
    young = 310.68E9;
    usf = 1.44339;
    An = usf;
  }
  else if (strcmp(style,"MoNbTi") ==  0){
    a=3.2373185E-10;
    b= sqrt(3.0)*a/2.0;
    mu = 32.30E9;//from ML potential //41.29E9; (DFT)
    young = 90.70E9;//ML potential //114.78E9; (DFT)
    usf = 0.761408228239499;
    An = usf;
  }
  else if (strcmp(style,"TaNbTi") ==  0){
    a=3.3014008E-10;//3.225E-10;
    b= sqrt(3.0)*a/2.0;
    mu = 29.55E9;//ML potential 
    young = 82.71E9;//ML potential
    usf = 0.549789143600001;
    An = usf;
  }
  else if (strcmp(style,"Nb_ML") == 0){
    a=3.326E-10;
    b= sqrt(3.0)*a/2.0;
    mu=24.89E9;
    young=71.52E9;
    usf=0.77774;
    An=usf;
  }
  else if (strcmp(style,"Mo_ML") == 0){
    a=3.173E-10;
    b= sqrt(3.0)*a/2.0;
    mu=117.80E9;
    young=305.67E9;
    usf=1.2321;
    An=usf;
  }
  else if (strcmp(style,"Ta_ML") == 0){
    a=3.322E-10;
    b= sqrt(3.0)*a/2.0;
    mu=48.79E9;
    young=134.60E9;
    usf=0.76308;
    An=usf;
  }

}

void Material::read_mpea_file(const char *filename, int N1, int N2, int N3, int NS)
{
  FILE * usfe_file = fopen(filename, "r");
  if (usfe_file == NULL) {
    fprintf(stderr, "Cannot open input file %s\n", filename);
    exit(1);
  }
  int n1,n2,n3,ns;
  fscanf(usfe_file, "%d %d %d %d", &n1, &n2, &n3, &ns);

  if(N1 != n1 || N2 != n2 || N3 != n3 || NS != ns){
    fprintf(stderr, "MPEA core energy file size does not match box %d %d %d %d %d %d %d %d\n",N1,n1,N2,n2,N3,n3,NS,ns);
    exit(1);
  }

  A_mpea = new double[N1*n2*n3*ns];

  int i, i_temp, j_temp, k_temp, ss_temp, index;
  double usfe_temp;
  for (i = 0; i < N1*n2*n3*ns; i++) {
    fscanf(usfe_file, "%d %d %d %d %lf", &i_temp, &j_temp, &k_temp, &ss_temp, &usfe_temp);
    index = i_temp*N2*N3 + j_temp*N3 + k_temp + ss_temp*N1*N2*N3;
    A_mpea[index] = usfe_temp*1E-3; //convert to J/m^2
  }
}
