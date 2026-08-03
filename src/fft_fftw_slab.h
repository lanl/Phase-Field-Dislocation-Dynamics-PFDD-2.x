/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
   ------------------------------------------------------------------------- */

#ifdef FFT_CLASS
FftStyle(fftw_slab,FFTW_Slab)

#else

#ifndef PFDD_FFT_FFTW_SLAB_H
#define PFDD_FFT_FFTW_SLAB_H

#include "pointers.h"
#include "mpi.h"
#include "fft.h"
#include "dfftw_mpi.h"

namespace PFDD_NS {

  class FFTW_Slab : public FFT {
  public:

    fftwnd_mpi_plan plan, iplan;
    fftw_complex *data_fftw, *work, *data_real, *data_real_2, *temp_data, *data_core, *data_strain, *data_strain_NoExt, *work_strain, *vstrain_complex, *elastic_energy_den_vs;

    FFTW_Slab(class PFDD_C *, int, char **);
    ~FFTW_Slab();

    void create_plan();
    void init();
    void setup();
    void rotate_stress(); // rotate stress from local to global
    void frec();
    void stiffness();
    void stiffness_CTB();
    void greens_function();
    void Bmatrix();
    void Fmatrix();

    void fft_send(void *, int, MPI_Datatype, int, int, MPI_Comm);
    void fft_recv(void *, int, MPI_Datatype, int, int, MPI_Comm ,MPI_Status *);

    void init_loop();
    void internal_energy();
    void internal_energy_vs();
    void update_order_parameter();
    void update_vs();
    void prepare_next_itr();
    void temp_stor_data();
    void reset_data();
    void stressfree_strain();
    void total_average_strain();
    void strain();
    void stress();
    void resolSS_Schmid();
    void resolSS_non_Schmid();
    //grad
    void gradient();
    void gradient_fcc();
    void core_energy_bcc_perfect();
    void core_energy_perfect();
    void core_energy_extended();
    void core_energy_sine();
    void core_energy_pyrII(); // Claire changed 8/1/18: was core_energy_isf_usf
    void core_energy_USFE_angle_tau();
    void core_energy_mpea();
    void core_energy_lup1D();
    void core_energy_lup2D();
    void core_energy_extended11();
    void core_energy_extended_test();
    void core_energy_extended_Ag();
    void core_energy_extended_Au();
    void core_energy_extended_Cu();
    void core_energy_extended_Cu_gbtest();
    void core_energy_extended_Ir();
    void core_energy_extended_Ni();
    void core_energy_extended_Pd();
    void core_energy_extended_Pt();
    void core_energy_extended_Rh();
    void core_energy_extended_2Px();
    void project_core_energy();
    void allocate();
    void prep_forward();
    void forward_mode1();
    void forward_mode2();
    void prep_backward();
    void backward_mode1();
    void backward_mode2();
    void initial_sxtal_mode1();
    void initial_sxtal_NoLoop();
    void initial_sxtal_mode2();
    void initial_sxtal_1LOrthoFCC();
    void initial_sxtal_1LInclinedFCC();
    void initial_sxtal_1LNonOrthoFCC();
    void initial_sxtal_1LOrthoBCC();
    void initial_sxtal_2LOrtho();
    void initial_sxtal_StraightOrtho();
    void initial_sxtal_StraightNonOrtho();
    void initial_sxtal_3D3ShcpBasal(); //claire added 07/28/18 to pair with app_3d_hcp_basal
    void initial_sxtal_3D3ShcpNotch(); //claire added 08/09/18 to pair with app_3d3s_hcp_notch
    void initial_sxtal_frank_read();
    void initial_sxtal_3SBCC();
    void initial_sxtal_edge();
    void initial_sxtal_loop();
    void initial_sxtal_cubic_fcc_edge();
    void initial_sxtal_cubic_fcc_screw();
    void initial_sxtal_3SFCC();
    void initial_sxtal_3SFCCPlanar();
    void initial_sxtal_6SFCCPlanar();
    void initial_sxtal_6SFCCPlanarKink();
    void initial_sxtal_6SFCCGPZone();
    void initial_sxtal_6SFCCPlanar_Obs();
    void initial_sxtal_6SFCCPlanar_ObsInf();
    void set_phase_2Px();
    void set_phase_3Px();

    void peierls(int);
    void screw_center();
    void xslip_sxtal_6SFCCPlanar();
    void xslip_check();
    void barrier_xslip();
    void energy_calc();

    void gb_energy_proportional();
    void gb_energy_extended();
    void no_gb_energy();
  };

}

#endif
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
