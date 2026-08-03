# Phase-Field-Dislocation-Dynamics-(PFDD)
Phase Field Dislocation Dynamics (PFDD) 

Project Scope:

Phase Field Dislocation Dynamics (PFDD) is a phase field approach for tracking the motion and interactions of individual dislocations in materials, mostly face-centered cubic (fcc) metals. The phase field approach relies on scalar order parameters or phase field variables to track a quantity of interest.  In PFDD, these phase field variables represent the passage of perfect dislocations across each active slip system, with one order parameter defined in each slip system.  Partial dislocations are represented with linear combinations of the order parameters.  The system is evolved through minimization of the total system energy, which is directly dependent on the phase field variables.  All of the physics is described by the total energy functional that is derived from a thermodynamic basis, and the system kinetics are evolved using a Ginzburg-Landau (GL) equation

Thus far, the codes here have been primarily used by a small group of researchers to study dislocation-dislocation interactions in metals, and the impact these interactions have on the overall deformation and material response of these materials.  Some related references including details about this formulation are:

1.	M. Koslowski, A.M. Cuitino, M. Ortiz, “A phase-field theory of dislocation dynamics, strain hardening and hysteresis in ductile single crystals”, Journal of the Mechanics and Physics of Solids, 50 (2002) 2597-2635.
2.	A. Hunter, I. J. Beyerlein, T. C. Germann, and M. Koslowski, “Influence of stacking fault energy surface on partial dislocations in fcc metals with a three-dimensional phase field dislocation dynamics model”, Physical Review B, 84 (2011) 144108.
3.	Y. Zeng, A. Hunter, I. J. Beyerlein, and M. Koslowski “A phase field dislocation dynamics model for a bicrystal interface system: An investigation into dislocation slip transmission across cube-on-cube interfaces”, International Journal of Plasticity, 79 (2016) 293-313.

Related review articles:

1.	A. Hunter, B. Leu, and I. J. Beyerlein, “A review of slip transfer: applications of mesoscale techniques”, Journal of Materials Science, 53 (2018) 5584-5603.
2.	I. J. Beyerlein and A. Hunter, “Understanding dislocation mechanics at the meso-scale using phase field dislocation dynamics”, Philosophical Transactions A, 374 (2016).

PFDD relies on a Fast Fourier Transform (FFT), which is built separately and not included in PFDD.  Any FFT solver can be employed, but as posted in the repo the parallel version is set-up to use [FFTW 2.1.5](https://www.fftw.org/download.html), and the serial version is set-up to use the Numerical Recipes fourn subroutine.  Note in order to run the codes, you need to obtain your own copies of either FFTW 2.1.5 or the Numerical Recipes fourn subroutine.

PFDD is used for investigating deformation in nanoscale (grain sizes of ~300 nm and less) materials, such as metals and alloys. This approach models the motion and interaction of individual defects, namely dislocations, in the material using scalar-valued phase field variables, also called order parameters. The system is evolved through energy minimization thus the model calculates the total energy density in terms of the phase field variables. The energy minimization is completed using the Ginzburg-Landau equation, and is implemented with explicit time integration. The total system energy can be comprised of several terms, including the strain energy (which describes dislocation-dislocation interactions), the energy due to an applied stress (dislocation interactions with the applied stress), and a core/lattice (perfect dislocations) or generalized stacking fault (partial dislocations) energy (described the dislocation core structure).  The latter term in particular may vary based on the crystal structure being modeled and is typically informed using lower length scale (e.g., atomistic) approaches, although no such (atomistic) calculations are completed within the PFDD algorithm. This basic formulation was previously reviewed by Los Alamos National Laboratory and released under license number C17113. This previously version PFDD v1.0 can be found on:
https://github.com/lanl/Phase-Field-Dislocation-Dynamics-PFDD

PFDD v2.0 is restructured and modularized from PFDD v1.0, and is now written in both C and C++ languages. New features include:
- new options of non-orthogonal computational grids;
- model Frank-Read sources;
- model dislocation cross-slip behavior;
- model multi-principal element alloys (MPEAs);
- non-Schmid effects in metals
- conductivity
- and new coupled formulation with micromechanical solvers

Some related references including details about new features are:

1. Peng, X., Hunter, A., Beyerlein, I. J., Lebensohn, R. A., Dayal, K., & Martinez, E. (2021). Non-orthogonal computational grids for studying dislocation motion in phase field approaches. Computational Materials Science, 200, 110834. https://doi.org/10.1016/j.commatsci.2021.110834

2. Kim, H., Mathew, N., Luscher, D. J., & Hunter, A. (2021). Phase field dislocation dynamics (PFDD) modeling of non-Schmid behavior in BCC metals informed by atomistic simulations. Journal of the Mechanics and Physics of Solids, 152, 104460. https://doi.org/10.1016/j.jmps.2021.104460

3. Fey, L. T., Hunter, A., & Beyerlein, I. J. (2022). Phase-field dislocation modeling of cross-slip. Journal of Materials Science, 1-15. https://doi.org/10.1007/s10853-021-06716-1

4. Blaschke, D. N., Miller, C., Mier, R., Osborn, C., Thomas, S. M., Tegtmeier, E. L., Winter, W. P., Carpenter, J. S., & Hunter, A. (2022). Predicting electrical conductivity in Cu/Nb composites: a combined model-experiment study. Journal of Applied Physics, 132, 045105. https://doi.org/10.1063/5.0096880 ([arXiv:2204.03777](https://www.arxiv.org/abs/2204.03777)).

5. Blaschke, D. N., Carpenter, J. S., & Hunter, A. (2024). Predicting electrical conductivity in bi-metal composites. Materials 17, 5049. https://doi.org/10.3390/ma17205049 ([arXiv:2409.04655](https://arxiv.org/abs/2409.04655)).

6. Brayan Murgas, Avanish Mishra, Nithin Mathew, Abigail Hunter; Phase field dislocation dynamics formulation coupled with Fourier based micromechanics solver and its application to grain boundary–dislocation interactions. J. Appl. Phys. 21 April 2026; 139 (15): 155101. https://doi.org/10.1063/5.0328415

PFDD relies on a Fast Fourier Transform (FFT), which is built separately and not included in PFDD. Any FFT solver can be employed, but as posted in the repo the current version is set-up to use [FFTW 2.1.5](https://www.fftw.org/download.html).

# Copyright
© 2022. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

# License
PFDD is distributed as open source software available under a [GPL3 license](GPLv3.pdf).

C22019

https://doi.org/10.11578/dc.20231117.5


***

# Compiling instructions Linux

Several Makefiles are distributed in the folder src/MAKE/. Depending on the user's system, minor modifications may be necessary.
For example, to compile on a linux machine using openmpi (after installing fftw 2.1.5, openmpi, and a C-shell such as csh or tcsh), run </br>
`cd src`
</br>
`make openmpi`
</br>
If fftw is not found or in an unusual location, set environment variable `FFTW` (such as `export FFTW = $(HOME)/fftw2` in bash)
</br>
The binary will be located in the src/ folder.

Example Debian based Linux distribution:
</br>
`apt install libopenmpi-dev fftw-dev fftw2 tcsh`
</br>
Note, the version of fftw2 shipped with Debian omits the 'd' prefix in the library names for double precision and the following edits are necessary if this version of fft2 is to be used:
</br>
in `src/fft_fftw_slab.h` and `src/fft_fftw_slab.cpp` replace `#include "dfftw_mpi.h"` with `#include "fftw_mpi.h"`
</br>
and in `src/MAKE/Makefile.openmpi` replace
</br>
`LIB = -ldfftw_mpi -ldfftw -lm -lstdc++`
</br>
with
</br>
`LIB = -lfftw_mpi -lfftw -lm -lstdc++`
</br>
then</br>
`cd src`
</br>
`make openmpi`

# Running

A few example input files are provided in the src folder, prefixed by "in."
</br>
To run an openmpi build on just one cpu with default settings on a linux machine call:</br>
`pfdd_openmpi -in <inputfile>`
</br>
Parallel runs must use n cpus where the number of regions defined in the inputfile divided by n is an integer, e.g.
</br>
`mpirun -n 4 pfdd_openmpi -in <inputfile>`
</br>
where the "region" command in `<inputfile>` has a multiple of 4 as its first entry.
The program writes a log file 'log.pfdd' by default unless specified otherwise via command line argument `-log <logfilename>`.

***

## Compiling instructions

Compiling instructions for Mac and HPC Linux system. There are two main steps, 1) installation of FFTW and 2) installation of PFDD code. 

As a recommendation in your home directory create a folder named *workspace*, e.g., by using `mkdir workspace`.

**This is optional but recommended:** Go to the *workspace* folder (`cd workspace`) and create a subfolder for your simulations named *simus* or *simulations*.

Note, if you type on the terminal `pwd` in the *workspace* folder you'll get something like this: `/Users/youruser/workspace`. 

Clone this repository in the *workspace* folder by using `git clone ssh://git@git.com/pfdd.git`. **The link `ssh://git@git.com/pfdd.git` is just an example**, you should copy the ssh link going to the repository, click on *code*, click on *code with SSH* and use the link provided there. For the github public version of the code this will be: `git@github.com:lanl/Phase-Field-Dislocation-Dynamics-PFDD-2.x.git`.

Note: In Linux based systems the ~/.ssh/config file does not require to have the `UseKeychain yes` which is a MacOS command. Also, you can add `IgnoreUnknown UseKeychain`.

---

### Installing PFDD on a HPC machine

#### Installing the fftw library

The following commands to load the modules are not necessary but you should check if they are available by running `module list`:

```
Currently Loaded Modules:
  1) craype-x86-spr           6) conduit/conduit-x86_64  11) intel/2023.2.0
  2) libfabric/1.22.0         7) gufi/gufi-x86_64        12) cray-mpich/8.1.28
  3) craype-network-ofi       8) cpe/23.12               13) cray-dsmml/0.2.2
  4) perftools-base/23.12.0   9) xpmem/2.9.7             14) cray-libsci/23.12.5
  5) slurm                   10) craype/2.7.30           15) PrgEnv-intel/8.5.0
```

The necessary modules are `intel` and `mpich`. Which are specified in the module 11 and 12 of the list above. These modules were loaded by default. The versions may change as they are updated, for instance, it may be `intel/2023.2.0` and `cray-mpich/8.1.28`.

In case the modules are not loaded, you can get the updated modules by running `module avail intel` or `module avail cray-mpich` to get the udpated version. Or `module spider intel` and `module spider cray-mpich` to see how to load the modules. By running `module spider intel` you'll get

```
----------------------------------------------------------------------------------
  intel:
----------------------------------------------------------------------------------
     Versions:
        intel/2023.1.0
        intel/2023.2.0
        intel/2024.1
     Other possible modules matches:
        PrgEnv-intel  intel-advisor  intel-ccl  intel-classic  intel-dal  ...

----------------------------------------------------------------------------------
  To find other possible module matches execute:

      $ module -r spider '.*intel.*'

----------------------------------------------------------------------------------
  For detailed information about a specific "intel" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider intel/2024.1
----------------------------------------------------------------------------------
```

running `module spider intel/2023.2.0` gives:

```
----------------------------------------------------------------------------------
  intel: intel/2023.2.0
----------------------------------------------------------------------------------

    You will need to load all module(s) on any one of the lines below before the "intel/2023.2.0" module is available to load.

      cpe/23.12  PrgEnv-intel/8.5.0
```

running `module spider cray-mpich/8.1.28` gives: 

```
----------------------------------------------------------------------------------
  cray-mpich: cray-mpich/8.1.28
----------------------------------------------------------------------------------

    You will need to load all module(s) on any one of the lines below before the "cray-mpich/8.1.28" module is available to load.

      cpe/23.12  PrgEnv-cray/8.5.0  cce/17.0.0
      cpe/23.12  PrgEnv-gnu/8.5.0  gcc/12.2.0
      cpe/23.12  PrgEnv-gnu/8.5.0  gcc/12.3.0
      cpe/23.12  PrgEnv-gnu/8.5.0  gcc/13.2.0
      cpe/23.12  PrgEnv-intel/8.5.0  intel-classic/2023.1.0
      cpe/23.12  PrgEnv-intel/8.5.0  intel-classic/2023.2.0
      cpe/23.12  PrgEnv-intel/8.5.0  intel-oneapi/2023.1.0
      cpe/23.12  PrgEnv-intel/8.5.0  intel-oneapi/2023.2.0
      cpe/23.12  PrgEnv-intel/8.5.0  intel-oneapi/2024.1
      cpe/23.12  PrgEnv-intel/8.5.0  intel-oneapi/2025.1
      cpe/23.12  PrgEnv-intel/8.5.0  intel-oneapi/2025.2.0-sandbox
      cpe/23.12  PrgEnv-intel/8.5.0  intel/2023.1.0
      cpe/23.12  PrgEnv-intel/8.5.0  intel/2023.2.0
      cpe/23.12  PrgEnv-intel/8.5.0  intel/2024.1
```

In the second to last line you can see that *cray-mpich/8.1.28* requires `cpe/23.12  PrgEnv-intel/8.5.0  intel/2023.2.0`, and the two first modules are the same necessary to load *intel/2023.2.0*. To install them you need to install intel and them cray mpich by running
```
module load cpe/23.12 
module load PrgEnv-intel/8.5.0
module load intel/2023.2.0
module load cray-mpich/8.1.28
```

And finally, check the modules were loaded by running `module list`:

```
Currently Loaded Modules:
  1) craype-x86-spr           6) conduit/conduit-x86_64  11) PrgEnv-intel/8.5.0
  2) libfabric/1.22.0         7) gufi/gufi-x86_64        12) intel/2023.2.0
  3) craype-network-ofi       8) cpe/23.12               13) cray-libsci/23.12.5
  4) perftools-base/23.12.0   9) craype/2.7.30           14) cray-mpich/8.1.28
  5) slurm                   10) cray-dsmml/0.2.2
```

Now go to the *workspace* directory:

`cd ~/workspace`

Copy the fftw tar file in pfdd/other or download it from the fftw website into the workspace folder:

`cp ~/workspace/pfdd/other/fftw-2.1.5.tar ~/workspace/`

Untar the fftw directory and configure and compile fftw: 

`tar -xvf fftw-2.1.5.tar`

The tar file can be downloaded from the fftw website but there is a copy of the file in `pfdd/other`.

`cd fftw-2.1.5`

`./configure --enable-mpi --enable-type-prefix --enable-static --prefix /user/youruser/workspace/fftw-2.1.5 `

The path `/user/youruser/workspace/fftw-2.1.5` should be properly changed to the path where the fftw-2.1.5 folder is.

`make install`

If you decided to install everything somewhere else, you can find the path to your current directory using the “pwd” command. If you need to debug your installation and re-install, you can type “make clean”. 

If all goes well, the fftw compilation should produce a “lib” and “include” directory. Record the paths to these directories. Now, we will shift to the PFDD compilation. 

#### Installing the pfdd code

First, we need to specify the location of the FFTW library. In the Make file you will see that the location is set using a variable called `FFTW`.

if you followed the steps and folder structure you should set the following path in the terminal:

```
export FFTW=~/workspace/fftw-2.1.5
```

You can check the Make file by running the following commands in th terminal:
`cd pfdd/MAKE`

`nano Makefile.openmpi`

locate the FFTW_INC and FFTW_LIB variables, the FFTW library you just defined is defined in the following lines of the Make file: 

```
PFDD_INC = -I/usr/local/include
PFDD_LIB = -L/usr/local/lib

FFTW_INC = -I/${FFTW}/include
FFTW_LIB = -L/${FFTW}/lib
```

The whole Makefile.openmpi is supposed to look like this

```
# linux = RedHat Linux box, Intel icc, MPICH2, FFTW

SHELL = /bin/sh

# ---------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =            icc
CC =            pgCC
CC =            mpicxx
#CCFLAGS =      -O2 \
                -funroll-loops -fstrict-aliasing -Wall -W -Wno-uninitialized
CCFLAGS =       -g
#CCFLAGS =      -O3
DEPFLAGS =      -M
LINK =          icc
LINK =          pgCC
LINK =          mpicxx
#LINKFLAGS =    -O2 \
                -funroll-loops -fstrict-aliasing -Wall -W -Wno-uninitialized
LINKFLAGS =     -g
#LINKFLAGS =    -O3
LIB =       -ldfftw_mpi -ldfftw -lm -lstdc++
ARCHIVE =       ar
ARFLAGS =       -rc
SIZE =          size

# ---------------------------------------------------------------------
# LAMMPS-specific settings
# specify settings for LAMMPS features you will use

# LAMMPS ifdef options, see doc/Section_start.html

PFDD_INC = -I/usr/local/include
PFDD_LIB = -L/usr/local/lib

FFTW_INC = -I/${FFTW}/include
FFTW_LIB = -L/${FFTW}/lib

# MPI library, can be src/STUBS dummy lib
# INC = path for mpi.h, MPI compiler settings
# PATH = path for MPI library
# LIB = name of MPI library

MPI_INC =
MPI_PATH =
MPI_LIB =

# JPEG library, only needed if -DLAMMPS_JPEG listed with LMP_INC
# INC = path for jpeglib.h
# PATH = path for JPEG library
# LIB = name of JPEG library

JPG_INC =
JPG_PATH =
#JPG_LIB =      /usr/local/lib/libjpeg.a
JPG_LIB =

# ---------------------------------------------------------------------
# build rules and dependencies
# no need to edit this section

EXTRA_INC = $(PFDD_INC) $(MPI_INC) $(JPG_INC) ${FFTW_INC}
EXTRA_PATH = $(MPI_PATH) $(JPG_PATH)
EXTRA_LIB = $(PFDD_LIB) $(MPI_LIB) $(JPG_LIB) ${FFTW_LIB}

# Link target

$(EXE): $(OBJ)
        $(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(EXTRA_LIB) $(LIB) -o $(EXE)
        $(SIZE) $(EXE)

# Library target

lib:    $(OBJ)
        $(ARCHIVE) $(ARFLAGS) $(EXE) $(OBJ)

# Compilation rules

%.o:%.cpp
        $(CC) $(CCFLAGS) $(EXTRA_INC) -c $<

%.d:%.cpp
        $(CC) $(CCFLAGS) $(EXTRA_INC) $(DEPFLAGS) $< > $@

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
include $(DEPENDS)
```

go back to the main pfdd folder 

`cd ..`

`make clean-all`

`make -j4 openmpi`

The -j4 flag just uses four cores to speed things up, you can also run “make openmpi”. If things go well, you should now see a pfdd_openmpi executable. Now we can try running an initial example. 

Copy a file in the simus folder. Create a running file `touch run.sh` and open it `nano run.sh`.

The structure of your batch submission script will look something like this: 
For Snow:

```
#!/bin/bash
#SBATCH -t 01:00:00     # walltime
#SBATCH -N 1           # number of cluster nodes
#SBATCH -n 4     # number of total tasks
#SBATCH --qos=debug #right now, this is set-up to run on debug queue so that it executes immediately
#SBATCH --reservation=debug
#SBATCH --partition=debug

module load intel/19.1.3 openmpi/3.1.6

mpirun -np 4 ~/workspace/pfdd/pfdd_openmpi -in in.fcc_plane1 > output/out-fcc-test
```

From a test around June 23rd 2025, this run file works:
```
#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=##
#SBATCH --time=12:00:00
cd $SLURM_SUBMIT_DIR
 
 
module load intel/19.1.3
module load cray-mpich/8.1.25
export LANG=en_US.UTF-8
 
ulimit -s unlimited
srun -n $SLURM_NTASKS ./pfdd_openmpi -in [INPUT FILE] > [OUTPUT FILE]
```
the additional modules are needed.

Or:

```
#!/bin/bash
#SBATCH -t 01:00:00     # walltime,
#SBATCH -N 1           # number of cluster nodes
#SBATCH -n 4     # number of total tasks
#SBATCH --qos=debug
#SBATCH --reservation=debug
#SBATCH --partition=debug

module load intel/19.1.3 cray-mpich/8.1.26 
#module load intel openmpi
cd .. 
srun -n 4 workspace/pfdd/pfdd_openmpi -in in.fcc_zerostress > output/out-fcc-zerostress-test
```

or using 

`srun -np $SLURM_NPROCS workspace/pfdd/pfdd_openmpi -in input`

Note the structure of the last line: 

`mpirun -np num_cores ./pfdd_executable_path -in input_script_path > output_file_path`

`sbatch run.sh`

***

---

### Installing PFDD locally on a Mac laptop 

#### Installing the fftw library

**Attention!** The following instructions have been tested on MacOS Monterey, Sonoma, Sequoia and Tahoe. They require admins rights on your computer to use `sudo`.

Open MPI is necessary for the installation. Please check if it is installed

```
mpirun --version
```

if it returns the following message

```
command not found: mpirun
```

please install open mpi

```
brew install open-mpi
```

and check again the version using `mpirun --version`, you should get a message that looks like this:

```
mpirun (Open MPI) 5.0.7
```

To install the fftw library copy the zip or tar file in the folder other inside the pfdd code and unzip/untar it in the *workspace* folder. This will create a folder named *fftw-2.1.5*. 

`cd fftw-2.1.5`

`sudo ./configure --enable-mpi --enable-type-prefix --enable-static --build=i686-apple-darwin17.7.0 --prefix /Users/youruser/workspace/fftw-2.1.5`

Note that the build flag above is necessary for compatibility of the old fftw version currently used in the code with building on a Mac with Apple silicon.

`sudo make install`

To successfully compile, you probably have to copy the lib and include files to /usr/local/ because they could not be successfully identified in the /Users/ directory. You can update the paths in my pfdd/MAKE/Makefile.openmpi file accordingly.

If all goes well, the fftw compilation should produce a “lib” and “include” directory inside the folder *~/workspace/fftw-2.1.5/*. Record the paths to these directories. 

#### Installing the pfdd code

Now, we will shift to the PFDD compilation. 

First, we need to specify the location of the FFTW library. In the Make file you will see that the location is set using a variable called `FFTW`.

if you followed the steps and folder structure you should set the following path in the terminal:

```
export FFTW=~/workspace/fftw-2.1.5
```

You can check the Make file by running the following commands in th terminal:
`cd pfdd/MAKE`

`nano Makefile.openmpi`

locate the FFTW_INC and FFTW_LIB variables, the *FFTW* variable you just defined is defined in the following lines of the Make file: 

```
PFDD_INC = -I/usr/local/include
PFDD_LIB = -L/usr/local/lib

FFTW_INC = -I/${FFTW}/include
FFTW_LIB = -L/${FFTW}/lib
```

Note that the paths `PFDD_INC`, `PFDD_LIB`, `FFTW_INC` and `FFTW_LIB` should point to the user's specific directory.

Compile the code

`cd pfdd`

`sudo make clean-all`

`sudo make openmpi` or `sudo make -j8 openmpi`

In the simulations file create another folder named test `mkdir test`, copy one of the test input files in the pfdd folder 

`cp ~/workspace/pfdd/input/inputfile ~/workspace/simus/test/`

`cd ~/workspace/simus/test`

Try to run the code:

`mpirun -np 1 ~/workspace/pfdd/pfdd_openmpi -in inputfile > out_sting_pfdd`

##### Debugging using lldb in MacOS:

You probably need to debug the code. You can use lldb to do so:

`lldb ~/workspace/pfdd/pfdd_openmpi`

setup a breakpoint

`b fft_fftw_slab.cpp:1000` 

run the simulation

`r -i inputfile`

the simulation will stop in your breakpoint or the next line with code. Use `n` or `c` to continue to the next line or the next time the breakpoint is executed. See the lldb documentation. You can use similar GDB in Linux for debugging.

---

## Note for new developers

 - Do not run simulations inside the pfdd folder that you cloned. The output files are too big and they are not necessary to track unless it is a unitary test (TODO).

 - Do not define material properties in other files but the material.cpp and material.h files. This will create issues when other people will use your code.

 - Do not modify the existing apps, fftw_* and solve_* files. If you are doing changes on these files it is better to create an additional file that fits your problem in the best way without interfering or modifying the existing code. 

 - Name variables, functions and files in a correct way.

 - Enjoy the code and science behind it!

## Input file structure:

### <span style="color:red">You can find more information in the manuals inside the documents folder.</span>.

An input script of PFDD looks like:

```
seed 576
app_style 1L1S100B core_energy 1 setup StraightOrtho
fft_style fftw_slab pbc 1 1 1 mode 1 primitive x 1 0 0 y 0 1 0 z 0 0 1
dimension 3
lattice sc/6n 1
region box block 0 64 0 64 0 64
create_box box
create_sites box
material preset Cu An 0.6
solve_style GL max_iter 40 tol 1e-4
diag_style strain
diag_style stress
temperature 700
scale 1.0 1.0 1.0
sigma initial 0.0 0.0 0.0 0.0 0.06 0.0 delta 0.0 0.0 0.0 0.0 0.0 0.0
stats 1.0e-8
dump 1 text 10 PFDD-Cu-rot-n64-S00.dump x y z xir1
run 1
```

This input script is for a simulation with a screw dislocation dipole in Cu in an orthogonal grid constructed with the basis x[100], y[010], z[001]. Detailed descriptions for each line in the input script are listed:

- **seed** seed for the random number generator. It is not required for the applications/examples presented in this manual but will be in the upcoming new applications.

- **app_style** 1L1S100B is the name of an application. Applications contain information to set up a dislocation network such as initial configuration, slip systems, as well as information required for a simulation such as time step and dislocation mobility. In this case, app_1L1S_100B.cpp assigns the initial configuration as a straight dislocation dipole along x axis with Burgers vector b = [100] and slip plane normal n = [001] (see method set_slip()). The initial configuration is specified by setup StraightOrtho. Time step and mobility are given in the declaration of class App1L1S100B. Some commonly used applications are included in the repository. Application usually has a default initial configuration but can be assigned with other preset initial configurations. See setup for details about different "application + setup" combinations.

- **core_energy** core energy flag. Currently available core energy forms are:
  1. perfect dislocation with sin2 approximation of the core energy, needs to specify An in material constants.
  5. line-character dependent core energy of BCC metals, needs to provide USF value of the chosen material.
  6. for non-Schmid application, same as the core energy ==1, except that the USFEmagnitude depends on the MRSSP angle and applied shear stress   .
  7. for BCC multi-principal element alloys. At each point the formulation is the same as the core energy == 1, except the USFE is variable depending on position. This core flag must be used with the mpea flag on the material line to read in an input script containing the local USFE values.

- **setup** would pick one preset initial configuration corresponding to the keyword. When there is no setup, application uses its default initial configuration (check setup_flag in application).\
Available configurations are listed with their keywords:
  - ***loop:*** a dislocation loop with radius 16b on yz-plane.
  - ***screw:*** a screw dislocation dipole on yz-plane assuming Burgers vector b is along the y axis; can be used with either orthogonal or non-orthogonal grid.
  - ***edge:*** an edge dislocation dipole on yz-plane specifically for BCC non-orthogonal grid, b has to be along the y axis. 
  - ***StraightOrtho:*** a dislocation dipole along y axis on xy-plane, depending on the slip system setup in the application, if b is along x - screw dipole; if b is along y - edge dipole. Can only be used for orthogonal grids.
  - ***StraightNonOrtho:*** an edge dislocation dipole on xy-plane specifically for FCC non-orthogonal grid, b has to be along x axis (for example application 1L1S100B).
  - ***frank_read:*** a Frank-Read source
  - ***none:*** no dislocation
  - ***mode1:*** 2 dislocation loops with radius 6b on 2 non-coplanar slip planes in the FCC orthogonal grid. When used with application 2L2S110B, the slip systems of these 2 loops are (111)[1 10] and (11 1)[1 10]. See the 2-loop example 3.1.3. 
  - ***mode2:*** 2 dislocation loops with radius 6b on 2 non-coplanar slip planes in the cubic grid. When used with application 2L2S110B, the slip systems of these 2 loops are (111)[1 10] and (11 1)[1 10].
  - ***1LOrthoFCC:*** a dislocation loop with radius 10b on xy-plane.
  - ***1LNonOrthoFCC:*** a dislocation loop with radius 10b on xy-plane and it is compatible with scaling of the primitive vectors (loop radius remains the same regardless of the scaling).
  - ***1LInclinedFCC:*** a dislocation loop with radius 10b on the (111) inclined slip plane of the cubic grid when used with application 1L1S110Binclined.
  - ***1LOrthoBCC:*** is the same as 1LOrthoFCC but compatible with the scaling of the primitive vectors.

##### Note that there are different ways to achieve the same dislocation setup. For example, 1L1S100B setup StraightOrtho with b along x axis gives a screw dislocation dipole on xy-plane, which can also be achieved if the line direction of StraightOrtho (ini101 tial_sxtal_StraightOrtho()) is along y and whichever application used has b along y. Therefore, there are flexible ways to make use of the above-mentioned applications and setups. Users are welcomed to explore other ways that are not covered here to suit their needs and preferences.

- **fft_style** PFDD calls for fft solver. While any FFT solver can be employed, the default one  106 is FFTW with its MPI library. FFTW MPI distribute data along the first dimension, i.e. the inputs and outputs are broken into disjoint blocks ("slabs") per process. 
  - ***fftw_slab*** stands for FFTW with MPI. 
  - ***pbc*** is followed with three integers to specify periodic boundary conditions along x; y; z directions, 1 - periodic and 0 - non-periodic. mode is followed with an integer to specify the way the data is organized for the fft calculation, 1 - mode 1 performs fft for each slip system separately and 2 - mode 2 performs fft all slip systems all at once. 
  - ***primitive*** corresponds to the primitive vectors on x; y; z directions, they are used to construct the computational grid. In this case, x 1 0 0 y 0 1 0 z 0 0 1 x[100], y[010], z[001]. Primitives for non-orthogonal grids can be found in the non-orthogonal examples. These are automatically normalized to unit vectors within the code.

- **dimension** PFDD code is set to model in a three-dimensional grid.

- **lattice** indicates the type of the lattice. For current applications, it remains sc/6n 1.

- **region** grid region along x; y; z directions. box block describe the region style. In 0 64 0 64 0 64, every two integers are a pair and stands for the start and the end along this direction, i.e. the grid goes from 0 to 64 on x; y; z.

- **create_box** and

- **create_sites** provide geometry information for building structured lattice with neighbors.

- **material** provides material constants for the simulation. It can either choose from the preset material list (see material.cpp) by having preset followed with the material name, or set material constants directly by having constant name with its value. Note that you can have both of them, but they can overwrite each other. To find out other preset materials and what constants you can assign, see material.cpp. In this case, it first calls material constants of Cu and then assign An to be 0.6. To simulate an MPEA, the mpea flag can be used followed by the name of a text file containing the local USFE values for each grid point. The first line should contain N1 N2 N3 NS, where N1, N2, and N3 are the dimensions of the cell and NS is the number of slip systems. Each subsequent line should read i j k s USFE, where USFE is the USFE on slip system s at point i, j, k.

- **solve_style** specifies the PFDD solver and stopping criteria of the simulation. GL - Ginzburg-Landau minimization. The only solver available now max_iter - max number of iteration in the minimization procedure. tol - tolerance for the convergence check. When max_iter and tol are both activated, simulations stops when it satisfies either of these two criteria.

- **diag_style** diagnostic style,
  - ***strain:*** includes average strain for every diagnostic output.
  - ***stress:*** includes average stress for every diagnostic output.
  - ***energy:*** includes average energy for every diagnostic output.

- **temperature** simulation temperature. It is not required for the applications/examples presented in this manual but will be in the upcoming new applications.

- **scale** is followed with three real numbers that are scaling factors for the three primitive vectors respectively. Primitive vectors are always normalized in the code upon assignment but can be scaled to desired length using scaling factors. In the case when they are 1.0 1.0 1.0, all primitives are of the unit length.

- **sigma** gives information about the applied stress.
  - ***initial:*** initial applied stress and
  - ***delta:*** increment of the applied stress
  The six numbers followed initial and delta are in the order of  xx;  yy;  zz;  xy;  xz;  yz. When using an orthogonal grid (either a rotated grid or a cubic grid), the applied stress follows the coordinate system as in "fft_style"; when using a non-orthogonal grid, the applied stress should be given on the cubic/Cartesian coordinate system. See Section 3.1 for more details.

- **stats** frequency for the diagnostic output.

- **dump** provides dump file information, including dump ID number, format, output interval, dump file name, names of the output items. Available output items include:
  - ***x,y,z*** the coordinates,
  - ***xir*** the real part of order parameters,
  - ***xii*** the imaginary part of order parameters,
  - ***pxx pyy pzz pxy pxz pyz*** stress components,
  - ***exx eyy ezz exy exz eyz*** strain components,
  - ***1 text 10 PFDD-Cu-rot-n64-S00.dump x y z xir1*** - the dump file 1, named PFDD-Cu-rot-n64-S00.dump, is a text file that contains output with an interval of every 10 iteration steps, and its output items are the coordinates x; y; z and the order parameter 1

## TODO:

 - Test cases.
 - Unitary tests.
 - CMake files.
