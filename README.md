# Phase Field Dislocation Dynamics (PFDD) version 2.x

PFDD is used for investigating deformation in nanoscale (grain sizes of ~300 nm and less) materials, such as metals and alloys. This approach models the motion and interaction of individual defects, namely dislocations, in the material using scalar-valued phase field variables, also called order parameters. The system is evolved through energy minimization thus the model calculates the total energy density in terms of the phase field variables. The energy minimization is completed using the Ginzburg-Landau equation, and is implemented with explicit time integration. The total system energy can be comprised of several terms, including the strain energy (which describes dislocation-dislocation interactions), the energy due to an applied stress (dislocation interactions with the applied stress), and a core/lattice (perfect dislocations) or generalized stacking fault (partial dislocations) energy (described the dislocation core structure).  The latter term in particular may vary based on the crystal structure being modeled and is typically informed using lower length scale (e.g., atomistic) approaches, although no such (atomistic) calculations are completed within the PFDD algorithm. This basic formulation was previously reviewed by Los Alamos National Laboratory and released under license number C17113. This previously version PFDD v1.0 can be found on:
https://github.com/lanl/Phase-Field-Dislocation-Dynamics-PFDD

PFDD v2.0 is restructured and modularized from PFDD v1.0, and is now written in both C and C++ languages. New features include:
new options of non-orthogonal computational grids;
model Frank-Read sources;
model dislocation cross-slip behavior;
model multi-principal element alloys (MPEAs);
and non-Schmid effects in metals.

Some related references including details about new features are:

1. Peng, X., Hunter, A., Beyerlein, I. J., Lebensohn, R. A., Dayal, K., & Martinez, E. (2021). Non-orthogonal computational grids for studying dislocation motion in phase field approaches. Computational Materials Science, 200, 110834. https://doi.org/10.1016/j.commatsci.2021.110834

2. Kim, H., Mathew, N., Luscher, D. J., & Hunter, A. (2021). Phase field dislocation dynamics (PFDD) modeling of non-Schmid behavior in BCC metals informed by atomistic simulations. Journal of the Mechanics and Physics of Solids, 152, 104460. https://doi.org/10.1016/j.jmps.2021.104460

3. Fey, L. T., Hunter, A., & Beyerlein, I. J. (2022). Phase-field dislocation modeling of cross-slip. Journal of Materials Science, 1-15. https://doi.org/10.1007/s10853-021-06716-1

4. Blaschke, D. N., Miller, C., Mier, R., Osborn, C., Thomas, S. M., Tegtmeier, E. L., Winter, W. P., Carpenter, J. S., & Hunter, A. (2022). Predicting electrical conductivity in Cu/Nb composites: a combined model-experiment study. Journal of Applied Physics, 132, 045105. https://doi.org/10.1063/5.0096880 ([arXiv:2204.03777](https://www.arxiv.org/abs/2204.03777)).

5. Blaschke, D. N., Carpenter, J. S., & Hunter, A. (2024). Predicting electrical conductivity in bi-metal composites. Materials 17, 5049. https://doi.org/10.3390/ma17205049 ([arXiv:2409.04655](https://arxiv.org/abs/2409.04655)).

PFDD relies on a Fast Fourier Transform (FFT), which is built separately and not included in PFDD. Any FFT solver can be employed, but as posted in the repo the current version is set-up to use [FFTW 2.1.5](https://www.fftw.org/download.html).

# Copyright
© 2022. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

# License
PFDD is distributed as open source software available under a [GPL3 license](GPLv3.pdf).

C22019

https://doi.org/10.11578/dc.20231117.5

# Compiling

Several Makefiles are distributed in the folder src/MAKE/. Depending on the user's system, minor modifications may be necessary.
For example, to compile on a linux machine using openmpi (after installing fftw 2.1.5, openmpi, and a C-shell such as csh or tcsh), run </br>
cd src
</br>
make openmpi
</br>
If fftw is not found or in an unusual location, edit src/MAKE/Makefile.openmpi and add a line pointing to it (such as FFTW = \$(HOME)/fftw2) above the line FFTW_INC = -I/\${FFTW}/include
</br>
The binary will be located in the src/ folder.

Example Debian based Linux distribution:
</br>
apt install libopenmpi-dev fftw-dev fftw2 csh
</br>
Note, the version of fftw2 shipped with Debian omits the 'd' prefix in the library names for double precision and the following edits are necessary if this version of fft2 is to be used:
</br>
in src/fft_fftw_slap.h and src/fft_fftw_slap.h.cpp replace #include "dfftw_mpi.h" with #include "fftw_mpi.h"
</br>
and in src/MAKE/Makefile.openmpi replace
</br>
LIB = -ldfftw_mpi -ldfftw -lm -lstdc++
</br>
with
</br>
LIB = -lfftw_mpi -ldfftw -lm -lstdc++
</br>
then</br>
cd src
</br>
make openmpi

# Running

A few example input files are provided in the src folder, prefixed by "in."
</br>
To run an openmpi build on just one cpu with default settings on a linux machine call:</br>
pfdd_openmpi -in <inputfile\>
</br>
Parallel runs must use n cpus where the number of regions defined in the inputfile divided by n is an integer, e.g.
</br>
mpirun -n 4 pfdd_openmpi -in <inputfile\>
</br>
where the "region" command in <inputfile> has a multiple of 4 as its first entry.
The programn's output is written to 'log.pfdd' by default unless specified otherwise via command line argument '-log <logfilename\>'.
