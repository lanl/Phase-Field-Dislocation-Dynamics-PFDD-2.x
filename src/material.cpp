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
#include "app.h"
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "material.h"
#include "fft.h"
#include "dfftw_mpi.h"
#include "memory.h"
#include "error.h"

#define WORD 256

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

Material::Material(PFDD_C *pfdd_p, int narg, char **arg, int N1, int N2, int N3, int NS, int NP) : Pointers(pfdd_p)
{
  // moved dslip up so gsfe lookup can use it
  // moved b into the preset material cases so it work for both BCC and FCC
  // b = sqrt(3.0)*(a)/2.0; //meters-->a is in meters
  // b = (a)/sqrt(2.0);
  dslip = 0.816; /*dslip=d/b and d=2(sqrt(3)a/2) or d = Constant*b*/
  omul = 1;
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
    else if (strcmp(arg[iarg],"lambda") == 0) {
      ll = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C11") == 0) {
      C11 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C12") == 0) {
      C12 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C44") == 0) {
      C44 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C13") == 0) {
      C13 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C14") == 0) {
      C14 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C15") == 0) {
      C15 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C16") == 0) {
      C16 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C22") == 0) {
      C22 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C23") == 0) {
      C23 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C24") == 0) {
      C24 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C25") == 0) {
      C25 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C26") == 0) {
      C26 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C33") == 0) {
      C33 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C34") == 0) {
      C34 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C35") == 0) {
      C35 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C36") == 0) {
      C36 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C45") == 0) {
      C45 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C46") == 0) {
      C46 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C55") == 0) {
      C55 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C56") == 0) {
      C56 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"C66") == 0) {
      C66 = atof(arg[iarg+1]);
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
    else if (strcmp(arg[iarg],"omul") == 0) {
      omul = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"Cn") == 0) {
      Cn = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"dlayer") == 0) {
      dlayer = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"Nlayers") == 0) {
      Nlayers = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"Ngrains") == 0) {
      Ngrains = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"Rgrain") == 0) {
      Rgrain = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"alpha_rho") == 0) {
      alpha_resistivity = atof(arg[iarg+1]);
      alpha_mfpel = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"pinterface") == 0) {
      pinterface = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"rho_disloc") == 0) {
      rho_disloc = atof(arg[iarg+1]); // TODO: will be removed in favor of automatic determination later
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"beta_disloc") == 0) {
      beta_disloc = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"El") == 0) {
      El[0] = atof(arg[iarg+1]); // set external electric field (required for conductivity calculations)
      El[1] = atof(arg[iarg+2]);
      El[2] = atof(arg[iarg+3]);
      iarg+=4;
    }
    else if (strcmp(arg[iarg],"preset") == 0) {
      set_constants(arg[iarg+1]);
      iarg+=2;
      strcpy(style,"preset");
    }
    else if (strcmp(arg[iarg],"mpea") == 0){
      read_mpea_file(arg[iarg+1], N1, N2, N3, NP);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"aflag") == 0) {   // read aniflag Aniso
      aflag = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"gsfe") ==0){
      read_gsfe_file(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"readHP") ==0){ // added this for reading activation energies from a file
      read_activation_file(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"error_tol") == 0) {  
      error_tol = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"rflag") == 0) {   // read reference flag 
      rflag = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"pbef") == 0) {   // read Peiersl barrier factor
      pbef = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"gbe") == 0) {   // read Grain boundary energy (inclusion)
      gbe = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"gbef") == 0) {   // read Peiersl barrier factor
      gbef = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"geometry") == 0) {   // read geometry flag 
      geometry = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"x0") == 0) {   // read origin flag 
      x0 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"y0") == 0) {   // read origin flag 
      y0 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"z0") == 0) {   // read origin flag 
      z0 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"xl") == 0) {   // read length flag 
      xl = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"yl") == 0) {   // read length flag 
      yl = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"zl") == 0) {   // read length flag 
      zl = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"sqp1") == 0) {   // read square parameter flag 
      sqp1 = atof(arg[iarg+1]);
      iarg+=2;
    }
    else if (strcmp(arg[iarg],"sqp2") == 0) {   // read square parameter flag 
      sqp2 = atof(arg[iarg+1]);
      iarg+=2;
    }
  }

//mu = 40.427E9;
//young = 110.50E9;
//ll = 111.20E9;

// this can be change in the input file adding 
// preset Cu aflag 1 ctol 5E-3

  // young = mu*(3*ll+2*mu)/(ll+mu);
  // xnu = young/2.0/(mu)-1.0;
  xnu = young/2.0/(mu)-1.0;
  ll = (C11+4*C12-2*C44)/5.0;
//  double kv;
//  kv = (C11+2.0*C12)/3.0;
//  ll = kv - 2.0*mu/3.0;


  /* virtual strain - aniso temp comment
  C44 = mu; //Pa
  nu = young/2.0/(mu)-1.0;
  C12 = 2.0*nu*C44/(1.0-2.0*nu); //Pa
  C11 = 2.0*C44+C12; //Pa
  */

  // moved b into the preset material cases so it work for both BCC and FCC
  // b = sqrt(3.0)*(a)/2.0; //meters-->a is in meters
  // b = (a)/sqrt(2.0);
  // dslip = 0.8165; /*dslip=d/b and d=2(sqrt(3)a/2) or d = Constant*b*/
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
  //  a = 3.52E-10;
  //  b = (a)/sqrt(2.0);
  //  // mu = 99.7E9;
  //  // young = 255.7E9;

  //  C11 = 246.5E9;
  //  C12 = 147.3E9;
  //  C44 = 127.4E9;
  //  mu = 96.28E9;
  //  young = 245.21E9;
  //  ll = 116.18E9;
   //mu = 49.6E9;
   //young = 136.3E9;
   // xnu = 0.374;
   // ll = 147.2E9;

   // Su et al. 2019
   a = 3.524E-10;
   b = (a)/sqrt(2.0);

   C11 = 247E9;
   C12 = 153E9;
   C44 = 122E9;
   mu = 92E9;
   young = 236.63E9;

   c0 = 503.144E-3;
   c1 = -106.533E-3;
   c2 = -83.741E-3;
   c3 = 28.777E-3;
   c4 = -2.916E-3;
   a1 = -123.917E-3;
   a3 = -38.533E-3;

  //  isf = 144.5E-3;  //25.0E-3;
  //  // usf = 289.0E-3; //270.0E-3;
  //  usf = 933.4E-3;
   
   condsigma = 1.43e7;
   mfpel = 5.87e-9; // mean free path of electrons
   alpha_resistivity = 0.006; // change in electrical resistivity per degree Kelvin close to 300K
   alpha_mfpel = 0.006; // change in mean free path with temperature
   Rgrain = 0.19; // model parameter: reflection of electrons at grain boundaries
   beta_disloc = 10.e-25; // electric resitivity per unit dislocation density [Ohm m^3]

  // Su et al. 2019
   isf = 152.13E-3;  //25.0E-3;
   usf = 300.89E-3; //270.0E-3;

   An = 801.47E-3;

 }
 else if (strcmp(style,"Al") == 0) {
  //  a = 4.04E-10;
   a = 4.05E-10; // Mishin 1999
   b = (a)/sqrt(2.0);

  //  mu = 28.4E9;
   mu = 26.26E9;
  //  young = 76.0E9;
   young = 75.87E9;
   ll = 59.269E9;

   c0 = 149.56346E-3;
   c1 = -4.23464E-3;
   c2 = -74.39331E-3;
   c3 = 33.84682E-3;
   c4 = -2.39438E-3;
   a1 = 67.75268E-3;
   a3 = -30.88995E-3;

   isf = 140.2E-3;
   usf = 177.0E-3;
   
   condsigma = 3.77e7;
   mfpel = 18.9e-9; // mean free path of electrons
   alpha_resistivity = 0.00429; // change in electrical resistivity per degree Kelvin close to 300K
   alpha_mfpel = 0.00429; // change in mean free path with temperature
   Rgrain = 0.17; // model parameter: reflection of electrons at grain boundaries
   beta_disloc = 1.5e-25; // electric resitivity per unit dislocation density [Ohm m^3]

  //  C11 = 114E9; // Mishin 1999
  //  C12 = 61.6E9; // Mishin 1999
  //  C44 = 31.6E9; // Mishin 1999

   C11 = 106.8E9;
   C12 = 60.4E9;
   C44 = 28.3E9;

   C22 = 106.8E9;
   C33 = 106.8E9;
   C55 = 28.3E9;
   C66 = 28.3E9;
   C13 = 60.4E9;
   C14 = 0.0;
   C15 = 0.0;
   C16 = 0.0;
   C23 = 60.4E9;
   C24 = 0.0;
   C25 = 0.0;
   C26 = 0.0;
   C34 = 0.0;
   C35 = 0.0;
   C36 = 0.0;
   C45 = 0.0;
   C46 = 0.0;
   C56 = 0.0;
 }
 else if (strcmp(style,"Au") == 0) {

  //  a = 4.17E-10;
  //  b = (a)/sqrt(2.0);
  //  mu = 18.8E9;
  //  young = 54.0E9;

  // Su et al. 2019
   a = 4.078E-10;
   b = (a)/sqrt(2.0);
   mu = 31.12E9;
   young = 88.04E9;
   C11 = 191E9;
   C12 = 162E9;
   C44 = 42.2E9;

   c0 = 161.832E-3;
   c1 = -40.659E-3;
   c2 = -17.320E-3;
   c3 = 3.197E-3;
   c4 = 0.472E-3;
   a1 = -64.718E-3;
   a3 = -13.282E-3;

  //  isf = 27.9E-3; //37.0E-3;
  //  usf = 66.5E-3;
   
   condsigma = 4.11e7;
   mfpel = 37.7e-9; // mean free path of electrons
   alpha_resistivity = 0.0034; // change in electrical resistivity per degree Kelvin close to 300K
   alpha_mfpel = 0.0034; // change in mean free path with temperature
   Rgrain = 0.17; // model parameter: reflection of electrons at grain boundaries
   beta_disloc = 2.6e-25; // electric resitivity per unit dislocation density [Ohm m^3]

  // Su et al. 2019
   isf = 26.79E-3; //37.0E-3;
   usf = 76.68E-3;
   An = 246.57E-3;
 }
 else if (strcmp(style,"Cu") == 0) {
   // a = 3.64E-10; //meters
   // b = (a)/sqrt(2.0);
   //mu = 54.5E9; //Pascal
  //  mu = 54.64E9; // From Janel's presentation 54.5
  //  young = 144.4E9; //Pascal
  //  ll = 100.64E9; // From Janel's presentation

  // Su et al. 2019
  a = 3.615E-10; // Experimental
  mu = 54.58E9; //Pascal
  young = 144.63E9; //Pascal
  C11=169E9;
  C12=122E9;
  C44=75.3E9;
   //
   // a = 3.61E-10; //Zeng bimetal (2016)
   // b = (a)/sqrt(2.0);
   // C11 = 168.4E9;
   // C12 = 121.4E9;
   // C44 = 75.4E9;
   // mu = 23.5E9; //Pascal
   // young = 66.7E9; //Pascal
   //
   // c0 = 344.948E-3;
   // c1 = -85.399E-3;
   // c2 = -40.421E-3;
   // c3 = 12.396E-3;
   // c4 = -0.674E-3;
   // a1 = -131.393E-3;
   // a3 = -19.816E-3;
   //
   // isf = 38.5E-3;
   // usf = 163.7E-3;
   
   condsigma = 5.95e7; // 5.96e7; // conductivity in SI units, redefine epsilonzero to use other units
   // conversion: 100% IACS (International Annealed Copper Standard) conductivity = 5.80e7 S/m (siemens per metre) in SI units
   alpha_resistivity = 0.00385; // change in electrical resistivity per degree Kelvin close to 300K
   mfpel = 39.9e-9; // mean free path of electrons
   alpha_mfpel = 0.00385; // change in mean free path with temperature
   Rgrain = 0.24; // model parameter: reflection of electrons at grain boundaries
   beta_disloc = 1.3e-25; // electric resitivity per unit dislocation density; some authors assume 2e-25 Ohm m^3

  //  a = 3.615E-10; // Mishin potential
   b = (a)/sqrt(2.0);
  //  mu = 55.2E9; // Mishin potential
  //  young = 138.4E9; // Mishin potential
  //  ll = 100.64E9; 
  //  C11 = 169.9E9;
  //  C12 = 122.6E9;
  //  C44 = 76.2E9;
   // Define the entire elasticity tensor for cases wher the crystals are rotated
  //  C22 = 169.9E9;
  //  C33 = 169.9E9;
  //  C55 = 76.2E9;
  //  C66 = 76.2E9;

  // Su et al. 2019
   C22 = 169E9;
   C33 = 169E9;
   C55 = 75.3E9;
   C66 = 75.3E9;
   C13 = 122E9;
   C14 = 0.0;
   C15 = 0.0;
   C16 = 0.0;
   C23 = 122E9;
   C24 = 0.0;
   C25 = 0.0;
   C26 = 0.0;
   C34 = 0.0;
   C35 = 0.0;
   C36 = 0.0;
   C45 = 0.0;
   C46 = 0.0;
   C56 = 0.0;

  //  c0 = 320.31E-3; // Mishin potential
  //  c1 = -73.07E-3;
  //  c2 = -46.28E-3;
  //  c3 = 13.36E-3;
  //  c4 = -0.27E-3;
  //  a1 = -109.72E-3;
  //  a3 = -22.90E-3;

   c0 = 344.75E-3;
   c1 = -85.42E-3;
   c2 = -40.36E-3;
   c3 = 12.46E-3;
   c4 = -0.70E-3;
   a1 = -131.28E-3;
   a3 = -19.70E-3;
  //  isf = 44.4E-3; // Mishin potential
  //  usf = 163.7E-3;

  // Su et al. 2019
   isf = 41.83E-3;
   usf = 160.52E-3;
   An = 519.05E-3;


 }
 else if (strcmp(style,"Pd") == 0) {
  //  a = 3.950E-10;
  //  b = (a)/sqrt(2.0);
  //  mu = 45.7E9;
  //  young = 125.8E9;

  // Su et al. 2019
   a = 3.89E-10;
   b = (a)/sqrt(2.0);
   mu = 52.48E9;
   young = 144.02E9;  
   C11=221E9;
   C12=171E9;
   C44=70.8E9;

   c0 = 308.65990E-3;
   c1 = -63.66363E-3;
   c2 = -48.76003E-3;
   c3 = 12.024463E-3;
   c4 = -1.11869E-3;
   a1 = -61.34161E-3;
   a3 = -23.33224E-3;

  //  isf = 138.1E-3;
  //  usf = 197.9E-3;

  // Su et al. 2019
   isf = 135.54E-3;
   usf = 214.93E-3;
   An = 476.98E-3;
 }
 else if (strcmp(style,"Ag") == 0) {
  //  a = 4.16E-10;
  //  b = (a)/sqrt(2.0);
  //  mu = 27.5E9;
  //  young = 74.8E9;

  // Su et al. 2019
   a = 4.086E-10;
   b = (a)/sqrt(2.0);
   mu = 33.3E9;
   young = 90.10E9;
   C11=122E9;
   C12=92E9;
   C44=45.5E9;

   c0 = 179.573E-3;
   c1 = -32.837E-3;
   c2 = -39.578E-3;
   c3 = 14.700E-3;
   c4 = -0.989E-3;
   a1 = -43.791E-3;
   a3 = -17.254E-3;

  //  isf = 17.8E-3; //37.0E-3;
  //  usf = 100.4E-3;
   
   condsigma = 6.3e7;
   mfpel = 53.3e-9; // mean free path of electrons
   alpha_resistivity = 0.0038; // change in electrical resistivity per degree Kelvin close to 300K
   alpha_mfpel = 0.0038; // change in mean free path with temperature
   Rgrain = 0.24; // model parameter: reflection of electrons at grain boundaries
   beta_disloc = 1.9e-25; // electric resitivity per unit dislocation density [Ohm m^3]

  // Su et al. 2019
   isf = 14.49E-3; //37.0E-3;
   usf = 92.01E-3;
   An = 307.39E-3;
 }
 else if (strcmp(style,"Pt") == 0) {
  //  a = 3.98E-10;
  //  b = (a)/sqrt(2.0);
  //  mu = 53.2E9;
  //  young = 149.0E9;

  // Su et al. 2019
   a = 3.923E-10;
   b = (a)/sqrt(2.0);
   mu = 65.1E9;
   young = 181.39E9;
   C11=347E9;
   C12=251E9;
   C44=76.5E9;

   c0 = 223.27233E-3;
   c1 = -15.38077E-3;
   c2 = -75.36747E-3;
   c3 = 21.70531E-3;
   c4 = -2.52876E-3;
   a1 = 58.60391E-3;
   a3 = -43.19019E-3;

  //  isf = 253.7E-3;
  //  usf = 257.5E-3;

  // Su et al. 2019
   isf = 280.87E-3;
   usf = 297.92E-3;
   An = 501.97E-3;
 }
 else if (strcmp(style,"Rh") == 0) {
  //  a = 3.84E-10;
  //  b = (a)/sqrt(2.0);
  //  mu = 154.1E9;
  //  young = 384.9E9;

  // Su et al. 2019
   a = 3.803E-10;
   b = (a)/sqrt(2.0);
   mu = 154.2E9;
   young = 387.92E9;
   C11=413E9;
   C12=194E9;
   C44=184E9;

   c0 = 556.35590E-3;
   c1 = -46.42878E-3;
   c2 = -183.45940E-3;
   c3 = 55.25737E-3;
   c4 = -5.04421E-3;
   a1 = 1.44561E-3;
   a3 = -72.22364E-3;

  //  isf = 184.6E-3;
  //  usf = 454.1E-3;

  // Su et al. 2019
   isf = 208.67E-3;
   usf = 475.34E-3;
   An = 958.3E-3;
 }
 else if (strcmp(style,"Ir") == 0) {
  //  a = 3.88E-10;
  //  b = (a)/sqrt(2.0);
  //  mu = 222.3E9;
  //  young = 549.7E9;

  // Su et al. 2019
   a = 3.839E-10;
   b = (a)/sqrt(2.0);
   mu = 221.2E9;
   young = 549.39E9;
   C11=580E9;
   C12=242E9;
   C44=256E9;

   c0 = 617.61574E-3;
   c1 = 3.28767E-3;
   c2 = -283.1158E-3;
   c3 = 96.12464E-3;
   c4 = -10.56298E-3;
   a1 = 157.20734E-3;
   a3 = -109.4680E-3;

  //  isf = 324.4E-3;
  //  usf = 614.9E-3;

  // Su et al. 2019
   isf = 385.87E-3;
   usf = 671.66E-3;
   An = 1250.63E-3;
 }
  else if (strcmp(style,"CoCrNi") == 0) { /* W.-R. Jian, Z. Xie, S. Xu, Y. Su, X. Yao, and I. J. Beyerlein, “Effects of lattice distortion and chemical short-range order on the mechanisms of deformation in medium entropy alloy CoCrNi,” Acta Materialia, vol. 199. Elsevier BV, pp. 352–369, Oct. 2020. doi: 10.1016/j.actamat.2020.08.044. */
   a = 3.556E-10;
   b = (a)/sqrt(2.0);
   mu = 71.84E9;
   young = 192.81E9;

   // NEEDS TO READ IN MPEA FILE FOR LATTICE ENERGY

   isf = -14.97E-3;
   usf = 310.34E-3;
   C11=252.6E9;
   C12=178.7E9;
   C44=95.1E9;
 }
 /* 11-term parameterizations for FCC metals, taken from Su et al., "Density functional theory calculations of generalized stacking fault energy surfaces for eight face-centered cubic transition metals", 2019.*/
 else if(strcmp(style,"Ag11") == 0) {
   a = 4.086E-10;
   b = (a)/sqrt(2.0);
   mu = 33.3E9;
   young = 90.10E9;

   c0 = 222.77E-3;
   c1 = -60.11E-3;
   c2 = -16.8E-3;
   c3 = 3.21E-3;
   c4 = 0.01E-3;
   c5 = -99.9E-3;
   c6 = -6.6E-3;
   c7 = 0.41E-3;
   c8 = 0.02E-3;
   c9 = -0.89E-3;
   c10 = -0.12E-3;

   isf = 14.49E-3; //37.0E-3;
   usf = 92.01E-3;
   C11=122E9;
   C12=92E9;
   C44=45.5E9;
   An = 307.39E-3;
 }
  else if(strcmp(style,"Au11") == 0) {
   a = 4.078E-10;
   b = (a)/sqrt(2.0);
   mu = 31.12E9;
   young = 88.04E9;

   c0 = 177.08E-3;
   c1 = -48.06E-3;
   c2 = -13.39E-3;
   c3 = 3.42E-3;
   c4 = -0.47E-3;
   c5 = -75.8E-3;
   c6 = -7.83E-3;
   c7 = 0.91E-3;
   c8 = 0.38E-3;
   c9 = -1.18E-3;
   c10 = -0.16E-3;

   isf = 26.79E-3; //37.0E-3;
   usf = 76.68E-3;
   C11 = 191E9;
   C12 = 162E9;
   C44 = 42.2E9;
   An = 246.57E-3;
 }
  else if(strcmp(style,"Cu11") == 0) {
   a = 3.615E-10; // Experimental
   // a = 3.64E-10; //meters
   b = (a)/sqrt(2.0);
   // mu = 54.5E9; //Pascal
   // young = 144.4E9; //Pascal
   mu = 54.58E9; 
   young = 144.63E9; 

   c0 = 383.24E-3;
   c1 = -107.47E-3;
   c2 = -22.86E-3;
   c3 = 3.82E-3;
   c4 = 0.07E-3;
   c5 = -172.98E-3;
   c6 = -9.5E-3;
   c7 = -0.04E-3;
   c8 = -0.12E-3;
   c9 = -0.95E-3;
   c10 = -0.27E-3;

   isf = 41.83E-3; // Mishin potential
   usf = 160.52E-3;
   C11=169E9;
   C12=122E9;
   C44=75.3E9;
   An = 519.05E-3;
 }
  else if(strcmp(style,"Ir11") == 0) {
   a = 3.839E-10;
   b = (a)/sqrt(2.0);
   mu = 221.2E9;
   young = 549.39E9;

   c0 = 932.83E-3;
   c1 = -185.45E-3;
   c2 = -131E-3;
   c3 = 11.95E-3;
   c4 = -0.77E-3;
   c5 = -189.38E-3;
   c6 = -33.78E-3;
   c7 = 5.95E-3;
   c8 = -0.53E-3;
   c9 = -9.72E-3;
   c10 = -0.6E-3;

   isf = 385.87E-3;
   usf = 671.66E-3;
   C11=580E9;
   C12=242E9;
   C44=256E9;
   An = 1250.63E-3;
 }
  else if(strcmp(style,"Ni-NM11") == 0) {
   a = 3.514E-10;
   b = (a)/sqrt(2.0);
   mu = 99.7E9;
   young = 220.1E9;

   c0 = 530.09E-3;
   c1 = -146.6E-3;
   c2 = -36.91E-3;
   c3 = 3.97E-3;
   c4 = -0.01E-3;
   c5 = -208.12E-3;
   c6 = -12.06E-3;
   c7 = 5.63E-3;
   c8 = -0.52E-3;
   c9 = -2.1E-3;
   c10 = -0.15E-3;

   isf = 139.01E-3;
   usf = 279.69E-3;
   C11 = 259.22E9;
   C12 = 172.74E9;
   C44 = 110.3E9;
 }
  else if(strcmp(style,"Ni11") == 0) {
   a = 3.524E-10;
   b = (a)/sqrt(2.0);
   mu = 92E9;
   young = 236.63E9;

   c0 = 589.67E-3;
   c1 = -161.16E-3;
   c2 = -42.21E-3;
   c3 = 5.07E-3;
   c4 = -0.18E-3;
   c5 = -229.7E-3;
   c6 = -16.07E-3;
   c7 = 4.29E-3;
   c8 = 0.03E-3;
   c9 = -2.08E-3;
   c10 = -0.12E-3;

   isf = 152.13E-3;  //25.0E-3;
   usf = 300.89E-3; //270.0E-3;

   C11 = 247E9;
   C12 = 153E9;
   C44 = 122E9;
   An = 801.47E-3;
 }
  else if(strcmp(style,"Pd11") == 0) {
   a = 3.89E-10;
   b = (a)/sqrt(2.0);
   mu = 52.48E9;
   young = 144.02E9;

   c0 = 360.73E-3;
   c1 = -93.69E-3;
   c2 = -27.91E-3;
   c3 = 1.74E-3;
   c4 = -0.14E-3;
   c5 = -117.38E-3;
   c6 = -9E-3;
   c7 = 2.82E-3;
   c8 = -0.45E-3;
   c9 = -2.31E-3;
   c10 = -0.18E-3;

   isf = 135.54E-3;
   usf = 214.93E-3;
   C11=221E9;
   C12=171E9;
   C44=70.8E9;
   An = 476.98E-3;
 }
  else if(strcmp(style,"Pt11") == 0) {
   a = 3.923E-10;
   b = (a)/sqrt(2.0);
   mu = 65.1E9;
   young = 181.39E9;

   c0 = 385.07E-3;
   c1 = -93.85E-3;
   c2 = -34.63E-3;
   c3 = 1.71E-3;
   c4 = -0.59E-3;
   c5 = -73.57E-3;
   c6 = -18.51E-3;
   c7 = 4.48E-3;
   c8 = -0.82E-3;
   c9 = -3.82E-3;
   c10 = -0.24E-3;

   isf = 280.87E-3;
   usf = 297.92E-3;
   C11=347E9;
   C12=251E9;
   C44=76.5E9;
   An = 501.97E-3;
 }
  else if(strcmp(style,"Rh11") == 0) {
   a = 3.803E-10;
   b = (a)/sqrt(2.0);
   mu = 154.2E9;
   young = 387.92E9;

   c0 = 708.48E-3;
   c1 = -151.06E-3;
   c2 = -90.14E-3;
   c3 = 9.2E-3;
   c4 = -0.3E-3;
   c5 = -189.92E-3;
   c6 = -22.71E-3;
   c7 = 3.04E-3;
   c8 = -0.43E-3;
   c9 = -5.38E-3;
   c10 = -0.79E-3;

   isf = 208.67E-3;
   usf = 475.34E-3;
   C11=413E9;
   C12=194E9;
   C44=184E9;
   An = 958.3E-3;
 }
 else if (strcmp(style,"Mg") == 0) {
   a = 3.347E-10;
   b = sqrt(3.0)*(a)/2.0; //BCC Mg
   // need to add Burges vector
   // b = * a
   mu = 42.456E9;
   ll = 37.66E9;
   young = 104.85E9;
   C11 = 84.68E9;
   C12 = 56.68E9;
   C44 = 61.4E9;

   c0 = 112.6205E-3;
   c1 = 0.7200E-3;
   c2 = -58.3781E-3;
   c3 = 28.1313E-3;
   c4 = -3.9212E-3;
   a1 = 25.2325E-3;
   a3 = -23.8545E-3;

   isf = 29.3216E-3;
   usf = 875.68E-3;

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
   
   condsigma = 7.634e6; // conductivity in SI units, redefine epsilonzero to use other units
   alpha_resistivity = 0.0035; // change in electrical resistivity per degree Kelvin close to 300K
   mfpel = 3.21e-9; // mean free path of electrons (Gall, JAP 2020)
   alpha_mfpel = 0.0035; // change in mean free path with temperature, estimate using CRC handbook
   Rgrain = 0.1; // model parameter: reflection of electrons at grain boundaries; TODO: find correct value (0.1 is a place holder to avoid nan)
   //~ beta_disloc = 0.0e-25; // electric resitivity per unit dislocation density [Ohm m^3]

   //  Estimate the Ta USFE (MRSSP angle, tau)--SNAP
   //fitted angle: -30^o to 30^o, stresses 0.25, 0.5, 0.75GPa
   a_slope = 51.59816407 ;
   a_b = -2.27173784 ;
   b_b = 0.9979452 ;
   c_slope = -4.16536 ;
   c_b = 1133.961364 ;
 }
 else if (strcmp(style,"W") == 0) {
   a = 3.1652E-10;
   b = sqrt(3.0)*(a)/2.0;
   mu = 160.6E9;
   young = 411.1E9;
   usf = 1850.E-3; // approx. value within range reported in https://doi.org/10.1016/j.commatsci.2021.110364
   
   condsigma = 1.855e7; // conductivity in SI units, redefine epsilonzero to use other units
   alpha_resistivity = 0.0042; // change in electrical resistivity per degree Kelvin close to 300K
   mfpel = 15.5e-9; // mean free path of electrons
   alpha_mfpel = 0.0042; // change in mean free path with temperature, estimate using CRC handbook
   Rgrain = 0.42; // model parameter: reflection of electrons at grain boundaries (Choi 2014)
   beta_disloc = 7.5e-25; // electric resitivity per unit dislocation density, taken from Brown 1977; Basinski, Dugdale, and Howie 1963 report 19e-25

 }
 else if (strcmp(style,"Cr") == 0) {
   a = 2.8848E-10;
   b = sqrt(3.0)*(a)/2.0;
   mu = 115.4E9;
   young = 279.3E9;
   //~ usf = 0.E-3;
   
   condsigma = 7.937e6; // conductivity in SI units, redefine epsilonzero to use other units
   alpha_resistivity = 0.003; // change in electrical resistivity per degree Kelvin close to 300K
   mfpel = 2.9e-9; // mean free path of electrons, taken from Idczak and Koscielniak Optica Applicata, Vol. X V I11, No. 2, 1988
   alpha_mfpel = 0.003; // change in mean free path with temperature, estimate using CRC handbook
   Rgrain = 0.1; // model parameter: reflection of electrons at grain boundaries; TODO: find correct value (0.1 is a place holder to avoid nan)
   //~ beta_disloc = 0.0e-25; // electric resitivity per unit dislocation density [Ohm m^3]

 }
 else if (strcmp(style,"Fe") == 0) {
  a = 2.8665E-10;
  b = sqrt(3.0)*(a)/2.0;
  mu = 87.2796E9;
  young = 225.1814E9;
  usf = 744.14E-3;

   condsigma = 1.041e6; // conductivity in SI units, redefine epsilonzero to use other units
   alpha_resistivity = 0.005; // change in electrical resistivity per degree Kelvin close to 300K
   mfpel = 40.0e-9; // mean free path of electrons
   alpha_mfpel = 0.005; // change in mean free path with temperature, estimate using CRC handbook
   Rgrain = 0.1; // model parameter: reflection of electrons at grain boundaries; TODO: find correct value (0.1 is a place holder to avoid nan)
   beta_disloc = 10.e-25; // electric resitivity per unit dislocation density, taken from Brown 1977;  Basinski, Dugdale, and Howie 1963 report 20e-25

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
    a = 3.347E-10;
    b = sqrt(3.0)*(a)/2.0;
    mu = 33.576E9;
    ll = 150.91E9;
    //mu = 43.89E9;
    young = 94.62E9;
    C11 = 249.01E9;
    C12 = 135.43E9;
    C44 = 18.1E9;
    usf = 676.775E-3;
    An = usf;
   
    condsigma = 6.334e6; // conductivity in SI units, redefine epsilonzero to use other units
    alpha_resistivity = 0.00318; // change in electrical resistivity per degree Kelvin close to 300K
    mfpel = 2.0e-9; // mean free path of electrons
    alpha_mfpel = 0.00318; // change in mean free path with temperature
    Rgrain = 0.45; // model parameter: reflection of electrons at grain boundaries
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
  // default values
  error_tol = 1E-6;
  rflag = 1; // Default value: Matrix as reference
  geometry = 0; // default geometry for inclusion: sphere
  x0 = y0 = z0 = 0;
  xl = yl = zl = 0;
  sqp1 = sqp2 = 2.0;
  pbef = 1.0; //Default Peiersl barrier factor for the inclusion
  gbe = 0.0; // Default Grain boundary energy
  gbef = 1.0; //Default Grain boundary energy factor for the inclusion
}
// Define new function to read in spatial parametrization of activation energy / prefactor along interface line
// example file in activation/1D_test.txt: first row: number of grid points along interface line, remaining rows: activation energy, prefactor
void Material::read_activation_file(const char *filename)
{

  FILE * activation_file = fopen(filename, "r");

  if (activation_file == NULL) {
    fprintf(stderr, "Cannot open input file %s\n", filename);
    exit(1);
  }

  int n1, i;
  double Hi, Pi;

  fscanf(activation_file, "%d", &n1);

  Hv = new double[n1];
  Hvg = new double[n1];
  Hvx = new double[n1];
  Pv = new double[n1];
  Pvg = new double[n1];
  Pvx = new double[n1];

  // Convention: read activation energies along interface line 
  for (i = 0; i < n1; i++) {
    fscanf(activation_file, "%lf %lf", &Hi, &Pi);
    Hvx[i] = Hi;
    Pvx[i] = Pi;
  }
  for (i = 0; i < n1; i++) {
    fscanf(activation_file, "%lf %lf", &Hi, &Pi);
    Hvg[i] = Hi;
    Pvg[i] = Pi;
  }

  fclose(activation_file);

}


void Material::read_mpea_file(const char *filename, int N1, int N2, int N3, int NP)
{
  // FILE * usfe_file = fopen(filename, "r");
  // if (usfe_file == NULL) {
  //   fprintf(stderr, "Cannot open input file %s\n", filename);
  //   exit(1);
  // }
  // int n1,n2,n3,ns;
  // fscanf(usfe_file, "%d %d %d %d", &n1, &n2, &n3, &ns);

  // if(N1 != n1 || N2 != n2 || N3 != n3 || NS != ns){
  //   fprintf(stderr, "MPEA core energy file size does not match box %d %d %d %d %d %d %d %d\n",N1,n1,N2,n2,N3,n3,NS,ns);
  //   exit(1);
  // }

  // A_mpea = new double[N1*n2*n3*ns];

  // int i, i_temp, j_temp, k_temp, ss_temp, index;
  // double usfe_temp;
  // for (i = 0; i < N1*n2*n3*ns; i++) {
  //   fscanf(usfe_file, "%d %d %d %d %lf", &i_temp, &j_temp, &k_temp, &ss_temp, &usfe_temp);
  //   index = i_temp*N2*N3 + j_temp*N3 + k_temp + ss_temp*N1*N2*N3;
  //   A_mpea[index] = usfe_temp*1E-3; //convert to J/m^2
  // }

  FILE * usfe_file = fopen(filename, "r");
  if (usfe_file == NULL) {
    fprintf(stderr, "Cannot open input file %s\n", filename);
    exit(1);
  }
  int n1,n2,n3,np;
  fscanf(usfe_file, "%d %d %d %d", &n1, &n2, &n3, &np);

  if(N1 != n1 || N2 != n2 || N3 != n3 || NP != np){
    fprintf(stderr, "MPEA core energy file size does not match box %d %d %d %d %d %d %d %d\n",N1,n1,N2,n2,N3,n3,NP,np);
    exit(1);
  }

  fcc7_mpea = new double[N1*n2*n3*np*7];

  int i, u_temp, i_temp, j_temp, k_temp, ss_temp, index;
  double sfe_temp;
  for (i = 0; i < 7*N1*n2*n3*np; i++) {
    fscanf(usfe_file, "%d %d %d %d %d %lf", &u_temp, &i_temp, &j_temp, &k_temp, &ss_temp, &sfe_temp);
    index = u_temp*N1*N2*N3 + i_temp*N2*N3 + j_temp*N3 + k_temp + ss_temp*7*N1*N2*N3;
    fcc7_mpea[index] = sfe_temp*1E-3; //convert to J/m^2
  }
}

void Material::read_gsfe_file(const char *filename)
{
  FILE * gsfe_file = fopen(filename, "r");
  if (gsfe_file == NULL) {
    fprintf(stderr, "Cannot open input file %s\n", filename);
    exit(1);
  }
  int i,j;
  if (app->core_flag == 8) { // 1D lookup table
    fscanf(gsfe_file, "%d", &i);
    Nx_lup = i;
    Ny_lup = 2;
  }
  else if (app->core_flag == 9) { // 2D lookup table
    fscanf(gsfe_file, "%d %d", &i, &j);
    Nx_lup = i;
    Ny_lup = j;
  }
  else {
    Nx_lup = 2;
    Ny_lup = 2;
  }
  if(fft->me == 0) printf("Nx_lup %d\nNy_lup %d\n",Nx_lup,Ny_lup);
  double gamma[Nx_lup*Ny_lup];

  int index;
  double xtemp, ytemp, gtemp;
  if(app->core_flag == 8){
    for (i = 0; i < Nx_lup; i++) {
      fscanf(gsfe_file, "%lf %lf", &xtemp, &gtemp);
      gamma[i] = gamma[Nx_lup + i] = gtemp * 1.E-3; // mJ to J
    }
  }
  else if(app->core_flag == 9) {
    for (j = 0; j < Ny_lup; j++) {
      for (i = 0; i < Nx_lup; i++) {
        fscanf(gsfe_file, "%lf %lf %lf", &xtemp, &ytemp, &gtemp);
        index = j*Nx_lup + i;
        gamma[index] = gtemp * 1.E-3; // mJ to J
      }
    }
  }
  else {
    fprintf(stderr, "GSFE lookup table is for GSFE lookup core energy functions only! Set core_flag to 8 (1D) or 9 (2D)");
    exit(1);
  }
  fclose(gsfe_file);
  compute_gsfe_derivatives(Nx_lup, Ny_lup, gamma);
}

void Material::compute_gsfe_derivatives(int Nx, int Ny, double *gsfe)
{
  int i, j, index;
  double sfe_i0[Ny], sfe_j0[Nx];
  double kx[Nx-1], ky[Ny-1], qx[Nx-1], qy[Ny-1];

  fftw_complex sfe[Nx-1][Ny-1], sfe_fft[Nx-1][Ny-1];
  fftw_complex sfe_fftqx[Nx-1][Ny-1], sfe_fftqy[Nx-1][Ny-1];

  fftw_real scale = (Nx-1)*(Ny-1);
  fftwnd_plan fwd_plan, bck_plan;

  fftw_complex sfe_fftdx[Nx-1][Ny-1], sfe_fftdy[Nx-1][Ny-1];
  sfe_dx = new double[(Nx-1)*(Ny-1)];
  sfe_dy = new double[(Nx-1)*(Ny-1)];
  for (j = 1; j < Ny; j++) {
    for (i = 1; i < Nx; i++) {
      sfe[i-1][j-1].re = 0.;
      sfe[i-1][j-1].im = 0.;
    }
  }

  for (i = 0; i < Nx; i++) {
    sfe_j0[i] = 0.;
  }

  for (j = 0; j < Ny; j++) {
    sfe_i0[i] = 0.;
  }

  for (j = 0; j < Ny; j++) {
    for (i = 0; i < Nx; i++) {
      index = j*Nx + i;
      if ((i != 0) && (j != 0)) {
        sfe[i-1][j-1].re = gsfe[index]/(mu*b*dslip);
      }
      else if (i == 0) {
        sfe_i0[j] = gsfe[index]/(mu*b*dslip);
      }
      else if (j == 0) {
        sfe_j0[i] = gsfe[index]/(mu*b*dslip);
      }
    }
  }
  for (j = 0; j < Ny-1; j++) {
    sfe[Nx-2][j].re = 0.5*(sfe[Nx-2][j].re + sfe_i0[j+1]);
  }

  for (i = 0;i < Nx-1; i++) {
    sfe[i][Ny-2].re = 0.5*(sfe[i][Ny-2].re + sfe_j0[i+1]);
  }

  double dx_lup, dy_lup;

  int Nx_lup1 = Nx - 1;
  int Ny_lup1 = Ny - 1;
  dx_lup = 1.0 / Nx_lup1;
  dy_lup = 1.0 / Ny_lup1;

  // initialize k and q
  double kx_lb, kx_ub, ky_lb, ky_ub;
  kx_lb = ceil((1 - Nx_lup1) / 2.0);
  kx_ub = floor(Nx_lup1 / 2.0);
  j = 0;
  for (i = 0; i < static_cast<int>(kx_ub) + 1; i++) {
    kx[j] = 2.0 * M_PI * i;
    j++;
  }
  for (i = static_cast<int>(kx_lb); i < 0; i++) {
    kx[j] = 2.0 * M_PI * i;
    j++;
  }
  if (j != Nx_lup1) {
    fprintf(stderr, "j is %d but should be %d\n", j, Nx_lup1);
    exit(1);
  }
  for (i = 0; i < Nx_lup1; i++) {
    qx[i] = sin(kx[i] * dx_lup) / dx_lup;
  }
  ky_lb = ceil((1 - Ny_lup1) / 2.0);
  ky_ub = floor(Ny_lup1 / 2.0);
  i = 0;
  for (j = 0; j < static_cast<int>(ky_ub) + 1; j++) {
    ky[i] = 2.0 * M_PI * j;
    i++;
  }
  for (j = static_cast<int>(ky_lb); j < 0; j++) {
    ky[i] = 2.0 * M_PI * j;
    i++;
  }
  if (i != Ny_lup1) {
    fprintf(stderr, "i is %d but should be %d\n", i, Ny_lup1);
    exit(1);
  }
  for (j = 0; j < Ny_lup1; j++) {
    qy[j] = sin(ky[j] * dy_lup) / dy_lup;
  }

  // Forward fft
  fwd_plan = fftw2d_create_plan(Nx_lup1, Ny_lup1, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
  fftwnd_one(fwd_plan, &sfe[0][0], &sfe_fft[0][0]);
  fftwnd_destroy_plan(fwd_plan);
  for (j = 0; j < Ny_lup1; j++) {
    for (i = 0; i < Nx_lup1; i++) {
      sfe_fftqx[i][j].re = -sfe_fft[i][j].im * qx[i] / scale;
      sfe_fftqx[i][j].im = sfe_fft[i][j].re * qx[i] / scale;
      sfe_fftqy[i][j].re = -sfe_fft[i][j].im * qy[j] / scale;
      sfe_fftqy[i][j].im = sfe_fft[i][j].re * qy[j] / scale;
  }
}

  // Backward fft
  bck_plan = fftw2d_create_plan(Nx_lup1, Ny_lup1, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
  fftwnd_one(bck_plan, &sfe_fftqx[0][0], &sfe_fftdx[0][0]);
  fftwnd_one(bck_plan, &sfe_fftqy[0][0], &sfe_fftdy[0][0]);
  fftwnd_destroy_plan(bck_plan);

  // Save real components as doubles
  for (j = 0; j < Ny_lup1; j++) {
    for (i = 0; i < Nx_lup1; i++) {
      index = j*Nx_lup1+i;
      sfe_dx[index] = sfe_fftdx[i][j].re;
      sfe_dy[index] = sfe_fftdy[i][j].re;
    }
  }
}
