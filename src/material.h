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

#ifndef PFDD_MATERIAL_H
#define PFDD_MATERIAL_H

#include "pointers.h"

namespace PFDD_NS {

  class Material : protected Pointers {
  public:
    char *style;                           // enum list of NONE,Al, Cu,etc

    double a;                    // lattice parameter
    double mu;                   // shear modulus
    double young;                // young modulus
    double xnu;                  // check???

    double c0,c1,c2,c3,c4,a1,a3; // parameters for the gamma surface
    double isf,usf;              // intrinsic and unstable stacking fault energies
    double An,Cn;
    double *A_mpea;            //variable An for an MPEA
    double aa0,aa1,aa2,aa3,aa4,bb1,bb2,bb3,bb4; // parameters for the hcp pyramidal2 gamma surface

    double C44,nu,C12,C11,ll;    // Elastic constants
    double S11,S12,S44;          // stiffness

    double b;                    // Burgers vector
    double dslip;                // interplanar distance

    double condsigma; // conductivity (TODO: make this a tensor)
    double alpha_resistivity = 0; // change in electrical resistivity per degree Kelvin above 300K (but close to it)
    double Rgrain=0.;
    double mfpel=1e-9; // mean free path of electrons
    double alpha_mfpel = 0; // change in mfp of electrons with temperature
    double dlayer=100e-9; // default average layer thickness
    int Nlayers = 1;
    int Ngrains = 1; // default number of grains per layer (in the transverse direction)
    double pinterface = 0.5; // default scattering probablility at an interface
    double epsilonzero=8.8541878128e-12; // vacuum permittivity in SI units according to CRC handbook
    double El[3]; // applied external electric field
    double beta_disloc = 0; // electric resitivity per unit dislocation density
    double rho_disloc = 0; // dislocation density, read from input file for now, 
    // but will eventually be automatically determined from the dislocations placed in the simulation box (TODO)
    
    double a_slope,a_b,b_b,c_slope,c_b; // paramters to determine USFE(X,tau)
    Material(class PFDD_C *, int, char **, int, int, int, int);
    ~Material();

  private:
    void set_constants(const char *);
    void read_mpea_file(const char *,int N1, int N2, int N3, int NS);
  };

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Lattice style does not match dimension

Self-explanatory.

E: Cannot use coloring without domain nx,ny,nz defined

UNDOCUMENTED

E: Color stencil is incommensurate with lattice size

Since coloring induces a pattern of colors, this pattern
must fit an integer number of times into a periodic lattice.

*/
