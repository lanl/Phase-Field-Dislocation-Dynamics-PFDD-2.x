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
    fftw_complex *data_fftw, *work, *data_real, *temp_data, *data_core, *data_strain, *work_strain;
    fftw_complex *conddata_fftw, *conddata_real, *conddata_bpart, *condtemp_data, *jcurrent, *jcurrent_fftw, *javerage;

    FFTW_Slab(class PFDD_C *, int, char **);
    ~FFTW_Slab();

    void create_plan();
    void init();
    void setup();
    void rotate_stress(); // rotate stress from local to global
    void frec();
    void stiffness();
    void greens_function();
    void Bmatrix();
    void Fmatrix();

    void fft_send(void *, int, MPI_Datatype, int, int, MPI_Comm);
    void fft_recv(void *, int, MPI_Datatype, int, int, MPI_Comm ,MPI_Status *);

    void init_loop();
    void internal_energy();
    void divcurrent();
    void update_order_parameter();
    void update_conduct_order_param();
    void prepare_next_itr();
    void prepare_next_conditr();
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
    void core_energy_bcc_perfect();
    void core_energy_perfect();
    void core_energy_extended();
    void core_energy_sine();
    void core_energy_pyrII(); // Claire changed 8/1/18: was core_energy_isf_usf
    void core_energy_USFE_angle_tau();
    void core_energy_mpea();
    void project_core_energy();
    void allocate();
    void prep_forward();
    void forward_mode1();
    void forward_mode2();
    void forwardcond();
    void prep_backward();
    void backward_mode1();
    void backward_mode2();
    void backwardcond();
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
