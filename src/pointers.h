/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
------------------------------------------------------------------------- */

// Pointers class contains ptrs to master copy of
//   fundamental LAMMPS class ptrs stored in lammps.h
// every LAMMPS class inherits from Pointers to access lammps.h ptrs
// these variables are auto-initialized by Pointer class constructor
// *& variables are really pointers to the pointers in lammps.h
// & enables them to be accessed directly in any class, e.g. atom->x

#ifndef PFDD_POINTERS_H
#define PFDD_POINTERS_H

#include "mpi.h"
#include "PFDD.h"
#include "Types.h"

namespace PFDD_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#define MAXLINE 1024
#define MAXWORD 128

  class Pointers {
  public:
  Pointers(PFDD_C *ptr) :
    pfdd_p(ptr),
      memory(ptr->memory),
      error(ptr->error),
      input(ptr->input),
      app(ptr->app),
      fft(ptr->fft),
      solve(ptr->solve),
      timer(ptr->timer),
      ranmaster(ptr->ranmaster),
      output(ptr->output),
      world(ptr->world),
      infile(ptr->infile),
      screen(ptr->screen),
      avgfile(ptr->avgfile),
      enerfile(ptr->enerfile),
      iterfile(ptr->iterfile),
      logfile(ptr->logfile) {}
    virtual ~Pointers() {}
    
  protected:
    PFDD_C *pfdd_p;
    Memory *&memory;
    Error *&error;
    Input *&input;
    App *&app;
    FFT *&fft;
    Solve *&solve;
    Timer *&timer;
    RanMars *&ranmaster;
    Output *&output;
    
    MPI_Comm &world;
    FILE *&infile;
    FILE *&screen;
    FILE *&logfile;
    FILE *&avgfile;
    FILE *&iterfile;
    FILE *&enerfile;
  };
  
}

#endif
