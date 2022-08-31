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

#ifndef PFDD_H
#define PFDD_H

#include "stdio.h"

namespace PFDD_NS {

  class PFDD_C {
  public:
    // ptrs to fundamental PFDD_C classes
    class Memory *memory;          // memory allocation functions
    class Error *error;            // error handling
    class Input *input;            // input script processing
    
    class App *app;                // application
    class FFT *fft;                // FFT solver
    class Solve *solve;            // update (GL/KMC/...)
    class Timer *timer;            // timer
    class RanMars *ranmaster;      // random number generator
    class Output *output;      // random number generator

    // ptrs to top-level LAMMPS-specific classes
    MPI_Comm world;                // MPI communicator
    int me;
    int nprocs;
    FILE *infile;                  // infile
    FILE *screen;                  // screen output
    FILE *logfile;                 // logfile
    
    char *suffix;                  // suffix to add to input script style names
    int suffix_enable;             // 1 if suffix enabled, 0 if disabled
    class Cuda *cuda;              // CUDA accelerator class
    
    PFDD_C(int, char **, MPI_Comm);
    ~PFDD_C();
    void create();
    void post_create();
    void destroy();
    
    void print_styles();
  };
  
}

#endif

/* ERROR/WARNING messages:

E: Invalid command-line argument

One or more command-line arguments is invalid.  Check the syntax of
the command you are using to launch LAMMPS.

E: Cannot use -reorder after -partition

Self-explanatory.  See doc page discussion of command-line switches.

E: Processor partitions are inconsistent

The total number of processors in all partitions must match the number
of processors LAMMPS is running on.

E: Must use -in switch with multiple partitions

A multi-partition simulation cannot read the input script from stdin.
The -in command-line option must be used to specify a file.

E: Can only use -pscreen with multiple partitions

Self-explanatory.  See doc page discussion of command-line switches.

E: Can only use -plog with multiple partitions

Self-explanatory.  See doc page discussion of command-line switches.

E: Cannot open universe screen file

For a multi-partition run, the master screen file cannot be opened.
Check that the directory you are running in allows for files to be
created.

E: Cannot open log.lammps

The default LAMMPS log file cannot be opened.  Check that the
directory you are running in allows for files to be created.

E: Cannot open universe log file

For a multi-partition run, the master log file cannot be opened.
Check that the directory you are running in allows for files to be
created.

E: Cannot open input script %s

Self-explanatory.

E: Cannot open screen file

The screen file specified as a command-line argument cannot be
opened.  Check that the directory you are running in allows for files
to be created.

E: Cannot open logfile

The LAMMPS log file named in a command-line argument cannot be opened.
Check that the path and name are correct.

E: Smallint setting in lmptype.h is invalid

It has to be the size of an integer.

E: Tagint setting in lmptype.h is invalid

Tagint must be as large or larger than smallint.

E: Bigint setting in lmptype.h is invalid

Size of bigint is less than size of tagint.

E: MPI_LMP_TAGINT and tagint in lmptype.h are not compatible

The size of the MPI datatype does not match the size of a tagint.

E: MPI_LMP_BIGINT and bigint in lmptype.h are not compatible

The size of the MPI datatype does not match the size of a bigint.

E: Small, tag, big integers are not sized correctly

See description of these 3 data types in src/lmptype.h.

E: Cannot use -cuda on without USER-CUDA installed

The USER-CUDA package must be installed via "make yes-user-cuda"
before LAMMPS is built.

*/
