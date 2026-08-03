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

#ifndef PFDD_FFT_H
#define PFDD_FFT_H

#include "pointers.h"
#include "material.h"
#include <complex.h>

namespace PFDD_NS {

  class FFT : protected Pointers {
  public:
    int me,nprocs;                    // proc info
    int procgrid[3];                  // assigned # of procs in each dim
    int user_procgrid[3];             // user request for procs in each dim
    int myloc[3];                     // which proc I am in each dim
    int procneigh[3][2];              // my 6 neighboring procs

    int box_exist;                    // 0 = not yet created, 1 = exists
    int dimension;                    // 1,2,3
    int nonperiodic;                  // 0 = periodic in all dims
    // 1 = non-periodic in any dim
    int xperiodic,yperiodic,zperiodic;  // 0 = non-periodic, 1 = periodic
    int periodicity[3];               // xyz periodicity as array

    double xprd,yprd,zprd;                               // global domain
    double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;    // global box bounds
    double subxlo,subxhi,subylo,subyhi,subzlo,subzhi;    // my portion of box
    int nx,ny,nz;               // mesh dimensions
    // only set for periodic sites from create_sites
    int mode;                   // Flag for the way the data is organized for the fft calculation
    int norder;                 // # of order parameters
    double **xi;             // order parameter for the phase field (xi)
    double **condxi;             // order parameter for the charge density phase field (condxi)
    double **xi_temp;
    double **xi_sum;             // order parameter for the phase field NOT CLEAR (xi_sum)
    double **dreg;             // Disregistry for the phase field (dreg)
    double **dreg_sum;             // Disregistry for the phase field NOT CLEAR (dreg_sum)
    double **dregp;             // Partial Disregistry for the phase field (dregp)
    double **dregp_sum;             // Partial Disregistry for the phase field NOT CLEAR (dregp_sum)
    double *xo;                 // order parameter for obstacles
    int *xcount;                // non-obstacle grid point count
    int eff_vol;                // non-obstacle effective volume
    double xinorm;              // norm of the change in xi
    int xflag;
    int peierlsz;
    int centerz1,centerz2,centerx1,centerx2;
    double peierlsdreg;
    double Hflag,X1flag,X2flag;
    int kmccheck1, kmccheck2;
    int **firstpass;
    double ximin;               //minimum order parameter in the slip plane i == 0
    double xiave;               // Average order parameter

    double intave, extave, latave;  // Average energies
    double obsden;              // obstacle density
    double *vstrain;            // order parameter for the virtual strain
    int *phase;                 // map of phases

    // Variable for LS micromechanics solver
    double errorconv;
    double convdot;
    
    int condnorder;           // # of order parameters for charge density
    double condxinorm;              // norm of the change in condxi
    double condxiave;               // Average charge density order parameter (DO WE NEED THIS?)
    double conductivity_global;  // calculated conductivity

    double *fx;                  // 1D array with all the frequencies in x
    double *fy;                  // 1D array with all the frequencies in y
    double *fz;                  // 1D array with all the frequencies in z
    double *f;                  // 1D array with all the frequencies NOT CLEAR
    double *r;                  // Not clear

    double ****C, ****S;               // stiffness matrix
    double ****C_CTB, ****S_CTB;       // stiffness matrix (CTB)
    double *****DC;             // Delta C for every phase 
    double *****Cinc;             // Delta C for every phase 
    double *****Cinvp;             // C inverse for every phase 
    double *****twoCC0invC0;             // C inverse for every phase 
    //double *****G;                 // green's function for each Fourier point
    double KCK[6][6], KCK_CTB[6][6];
    double RCR[3][3][3][3];     // Rotated 4 order stiffness matrix
    double *AA;
    double *BB;                 // Green's function
    double *FF, *PD, *PD2, *PDHet;                 // Green's function
    double *DD;                 // Green's function
    double *****MM;                 // virtual strain's
    double *EE;                 // virtual strain's
    // heterogeneous form
    double *AAhet;
    double *BBhet;              // Green's function
    double *DDhet;              // Green's function
    double *FFhet;                 // Green's function
    double *EEhet;                 // virtual strain's
    double *MGreen;             // Modified Green function multiplied by wave vector
    // heterogeneous form with second term
    double *AAhet2;
    double *BBhet2;              // Green's function
    double *DDhet2;              // Green's function
    double *FFhet2;                 // Green's function
    double *EEhet2;                 // virtual strain's
    //Grad
    double *gradx;              // gradient on x of xi
    double *grady;              // gradient on y of xi
    double *gradz;              // gradient on z of xi
    double *theta;              // line character angle
    double *hess;              // Hessian of xi
    double *nhn;                // 2nd order directional derivative of xi on grad direction
    double *picked;             // grid points picked for locating dislocation

    //Grad FCC
    double *gradx_d;              // gradient on x of xi
    double *grady_d;              // gradient on y of xi
    double *gradz_d;              // gradient on z of xi
    double *theta_d;              // line character angle
    double *hess_d;              // Hessian of xi
    double *nhn_d;                // 2nd order directional derivative of xi on grad direction
    double *picked_d;             // grid points picked for locating dislocation
    double *dt_kmc;

    double *fcore;              // 1D array for the gamma surface
    double *f_core;              // 1D array for the gamma surface. Added for MPI in extended core energy by CLaire 7/31/18
    double *df1core;            // derivative of the gamma surface in 1
    double *df2core;            // derivative of the gamma surface in 2
    double *df3core;            // derivative of the gamma surface in 3
    double *dE_core;            // sum of the derivatives of the gamma surface
    double E_core;              // core energy
    double *dE_gb;            // sum of the derivatives of the grain boudnary energy

    double *E_int;
    double *E_ext;
    double *E_lat;
    double *datareal;

    double conductivity;        // TODO: make this a tensor

    double *delta;              // 1D array for order paramter (projected)
    double *ddelta;             // 1D array for gradient of delta along x

    double **xn;                // glide plane normals
    double **xb;                // Burgers vectors

    double **sigma;             // Input stress
    double **deltasig;          // increment of sigma
    double **sigma_rot;         // stress in global coords
    double *data_sigma;         // output with the stress distribution
    double *data_asigma;         // output with the stress distribution
    double *data_tsigma;         // output with the stress distribution
    double *data_psigma;         // output with the stress distribution
    double *data_usigma;         // output with the stress distribution
    double *data_rsigma;         // output with the stress distribution
    double *data_sigma_NoExt;
    double *tau;                // Peach-Koehler force
    double *tau_2;
    double **ave_sigma, **ave_sigma_NoExt;          // average stress global
    double **local_sigma, **local_sigma_NoExt;       // average stress local

    double *Rij, *Fij, *Fijinv;                // Coordinate transformation
    double *Gamma;

    double ***eps;              // strain
    double *data_eps, *data_eps_NoExt;          // auxiliary strain
    double *data_epsd;          // auxiliary strain
    double **aveps, **aveps_NoExt;            // local average strain
    double **avepsd;            // average
    double **avepst;            // average
    double **avepsts;           // average
    double **ave_epsd;          // global average
    double **ave_eps, **ave_eps_NoExt;          // global average strain
    double **ave_Vs;          // global average
    double **aveVs;          // global average strain

    double ****avec;            // average
    double ****ave_c;            // average

    double *data_ieps;          // iterative strain

    int slip_systems;            // # of slip systems (NS)
    int stress_inc;              // # of stress increments (NSI)
    int num_planes;              // # of glide planes (NP)

    int local_x;                // # of local rows for the FFT
    int local_x_start;          // id of starting local row
    int local_ny_after_trans;      // if fftwnd_mpi is called with FFTW_TRANSPOSED_ORDER output
    // then y will be the first dimension of the output,
    // and the local y extent will be given by local_ny_after_trans
    int local_y_start_after_trans; // id of starting y
    int total_local_size;          // the number of fftw_complex elements that
    // you must allocate for your local data
    int local_y;                // # of local grid points in y
    int local_z;                // # of local grid points in z

    class Lattice *lattice;                  // user-defined lattice
    //~ class Material *material[5] {pfdd_p, pfdd_p, pfdd_p, pfdd_p, pfdd_p};            // user-defined materials (max: 5)
    class Material *material;            // user-defined material, always points to last one defined
    class Material *materialA;
    class Material *materialB;
    class Material *materialC;
    class Material *materialD;
    class Material *materialE;
    int Nmaterials=0;        // number of materials read from file

    int nregion;                             // # of defined Regions
    int maxregion;                           // max # list can hold
    class Region **regions;                  // list of defined Regions

    // New material and region input
    class Material *material2;                  // user-defined material
    class Material *material3;                  // user-defined material
    class Material *dislocation;                  // user-defined material

    int primitive;              // Flag 1 if lattice is primitive
    double vol;               // volume of the lattice vectors det(sclprim)
    double norvol;            // volume of the normalized lattice vectors det(prim)
    double prim[3][3];          // Matrix with the normalized primitive vectors
    double invprim[3][3];       //inverse of normalized primitive vectors
    double sclprim[3][3];       //scaled primitive vectors
    double four[3][3];          // Matrix with the scaled primitive vectors in Fourier space
    double normfour[3];         // Norm of the reciprocal vectors
    double xprdf,yprdf,zprdf;                               // global domain in rotated axis
    double boxxlof,boxxhif,boxylof,boxyhif,boxzlof,boxzhif;    // global box bounds
    double cntf[3];                    // Center point in rotated space

    double ortrot[3][3];          // orthonormal rotation matrix
    double invortrot[3][3];       //inverse of orthonormal rotation matrix
    double norvolort;

    FFT(class PFDD_C *, int, char **);
    ~FFT();

    void set_box();
    void set_lattice(int, char **);
    void set_material(int, char **);
    void add_region(int, char **);
    int find_region(char *);
    void set_boundary(int, char **);

    // New material and region input
    void set_material2(int, char **);
    void set_material3(int, char **);
    void set_dislocation(int, char **);

    void procs2domain_1d();
    void procs2domain_2d();
    void procs2domain_3d();

    double determinant(double [6][6], int);
    void rotate_stiffness(double [6][6], double [3], double [3], double [3]);
    void rotate_stiffness_CTB(double [6][6], double [3], double [3], double [3]);
    void rotate_stiffnessfour(double [3][3][3][3], double [3][3]);


    virtual void create_plan() = 0;
    virtual void init() = 0;
    virtual void setup() = 0;
    virtual void frec() = 0;
    virtual void stiffness() = 0;
    virtual void stiffness_CTB() = 0;
    virtual void greens_function() = 0;
    virtual void internal_energy() = 0;
    virtual void internal_energy_vs() = 0;
    //~ virtual void divcurrent() = 0;
    virtual void update_order_parameter() = 0;
    //~ virtual void update_conduct_order_param() = 0;
    virtual void update_vs() = 0;
    virtual void prepare_next_itr() = 0;
    //~ virtual void prepare_next_conditr() = 0;
    virtual void temp_stor_data() = 0;
    virtual void reset_data() = 0;
    virtual void Bmatrix() = 0;
    virtual void Fmatrix() = 0;
    virtual void allocate() = 0;
    virtual void prep_forward() = 0;
    virtual void forward_mode1() = 0;
    virtual void forward_mode2() = 0;
    //~ virtual void forwardcond() = 0;
    virtual void prep_backward() = 0;
    virtual void backward_mode1() = 0;
    virtual void backward_mode2() = 0;
    //~ virtual void backwardcond() = 0;
    virtual void initial_sxtal_mode1() = 0;
    virtual void initial_sxtal_mode2() = 0;
    virtual void initial_sxtal_NoLoop() = 0;
    virtual void initial_sxtal_2LOrtho() = 0;
    virtual void initial_sxtal_1LOrthoFCC() = 0;
    virtual void initial_sxtal_1LNonOrthoFCC() = 0;
    virtual void initial_sxtal_1LInclinedFCC()=0;
    virtual void initial_sxtal_1LOrthoBCC() = 0;
    virtual void initial_sxtal_StraightOrtho() = 0;
    virtual void initial_sxtal_StraightNonOrtho() = 0;
    virtual void initial_sxtal_3D3ShcpBasal() = 0; //claire added 07/28/18 to pair with app_3d_hcp_basal
    virtual void initial_sxtal_3D3ShcpNotch() = 0; //claire added 08/09/18 to pair with app_3d_3s_ortho_notch
    virtual void initial_sxtal_frank_read() = 0;
    virtual void initial_sxtal_3SBCC() = 0;
    virtual void initial_sxtal_edge() = 0;
    virtual void initial_sxtal_loop() = 0;
    virtual void initial_sxtal_cubic_fcc_edge() = 0;
    virtual void initial_sxtal_cubic_fcc_screw() = 0;
    virtual void initial_sxtal_3SFCC() = 0;
    virtual void initial_sxtal_3SFCCPlanar() = 0;
    virtual void initial_sxtal_6SFCCPlanar() = 0;
    virtual void initial_sxtal_6SFCCPlanarKink() = 0;
    virtual void initial_sxtal_6SFCCGPZone() = 0;
    virtual void initial_sxtal_6SFCCPlanar_Obs() = 0;
    virtual void initial_sxtal_6SFCCPlanar_ObsInf() = 0;
    virtual void set_phase_2Px() = 0;
    virtual void set_phase_3Px() = 0;

    virtual void fft_send(void *, int, MPI_Datatype, int, int, MPI_Comm) = 0;
    virtual void fft_recv(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *) = 0;


    virtual void gradient()= 0; //Grad
    virtual void gradient_fcc()= 0; //Grad
    virtual void core_energy_bcc_perfect() = 0;
    virtual void core_energy_perfect() = 0;
    virtual void core_energy_extended() = 0;
    virtual void core_energy_sine() = 0;
    virtual void core_energy_pyrII() = 0;//Claire changed 8/1/18: was core_energy_isf_usf
    virtual void core_energy_USFE_angle_tau() =0; //added HJ
    virtual void core_energy_mpea() = 0;
    virtual void core_energy_lup1D() = 0;
    virtual void core_energy_lup2D() = 0;
    virtual void core_energy_extended11() = 0;
    virtual void core_energy_extended_test() = 0;
    virtual void core_energy_extended_Ag() = 0;
    virtual void core_energy_extended_Au() = 0;
    virtual void core_energy_extended_Cu() = 0;
    virtual void core_energy_extended_Cu_gbtest() = 0;
    virtual void core_energy_extended_Ir() = 0;
    virtual void core_energy_extended_Ni() = 0;
    virtual void core_energy_extended_Pd() = 0;
    virtual void core_energy_extended_Pt() = 0;
    virtual void core_energy_extended_Rh() = 0;
    virtual void core_energy_extended_2Px() = 0;
    virtual void project_core_energy() = 0;

    virtual void init_loop() = 0;
    virtual void stressfree_strain() = 0;
    virtual void total_average_strain() = 0;
    virtual void strain() = 0;
    virtual void stress() = 0;
    virtual void resolSS_Schmid() = 0;
    virtual void resolSS_non_Schmid() = 0;
    virtual void rotate_stress() = 0;

    virtual void peierls(int) = 0;
    virtual void screw_center() = 0;
    virtual void xslip_sxtal_6SFCCPlanar() = 0;
    virtual void xslip_check() = 0;
    virtual void barrier_xslip() = 0;
    virtual void energy_calc() = 0;

    virtual void gb_energy_proportional() = 0;
    virtual void gb_energy_extended() = 0;
    virtual void no_gb_energy() = 0;

  };

}

#endif

/* ERROR/WARNING messages:

E: Box bounds are invalid

Lo bound >= hi bound.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Reuse of region ID

Self-explanatory.

E: Invalid region style

Self-explanatory.

E: Boundary command currently only supported by on-lattice apps

UNDOCUMENTED

E: App style proc count is not valid for 1d simulation

There can only be 1 proc in y and z dimensions for 1d models.

E: App style proc count is not valid for 2d simulation

There can only be 1 proc in the z dimension for 2d models.

*/
