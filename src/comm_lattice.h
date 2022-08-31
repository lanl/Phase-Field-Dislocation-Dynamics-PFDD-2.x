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

#ifndef PFDD_COMM_LATTICE_H
#define PFDD_COMM_LATTICE_H

#include "mpi.h"
#include "pointers.h"

namespace PFDD_NS {

class CommLattice : protected Pointers {
 public:
  CommLattice(class PFDD_C *); ~CommLattice();
  void init(int, int *);
  void all();
  void all_reverse();
  //void sector(int);
  //void reverse_sector(int);

 private:
  int me,nprocs;
  int delghost,delreverse;

  struct Swap {
    int nsend,nrecv;               // number of messages to send/recv
    int *sproc;                    // proc for each send message
    int *scount;                   // size of each send message in sites
    int *smax;                     // max size of each send message in sites
    int **sindex;                  // list of my lattice indices for each send
    int *sibuf;                    // biggest int send message
    double *sdbuf;                 // biggest double send message
    int *rproc;                    // proc for each recv message
    int *rcount;                   // size of each recv message in sites
    int *rmax;                     // max size of each recv message in sites
    int **rindex;                  // list of my lattice indices for each recv
    int **ribuf;                   // each int recv message
    double **rdbuf;                // each double recv message
    MPI_Request *request;          // MPI datums for each recv message
    MPI_Status *status;
  };

  struct Site {
    tagint id_global;
    int index_local;
    int proc;
  };

  Swap *allswap;
  Swap *reverseswap;
  //Swap **sectorswap;
  //Swap **sectorreverseswap;
  //int nsector;

  int ninteger,ndouble;
  int site_only;                     // only 1 int, no doubles
  int **iarray;
  double **darray;
  int *site;                         // simply points to iarray[0]

  Swap *create_swap_all();
  Swap *create_swap_all_reverse();
  //Swap *create_swap_sector(int, int *);
  //Swap *create_swap_sector_reverse(int, int *);
  void free_swap(Swap *);

  void create_send_from_list(int, Site *, Swap *);
  void create_send_from_recv(int, int, Site *, Swap *);
  void create_recv_from_send(int, int, Site *, Swap *);
  void create_recv_from_list(int, Site *, Swap *);

  void perform_swap_site(Swap *);
  void perform_swap_int(Swap *);
  void perform_swap_double(Swap *);
  void perform_swap_general(Swap *);
};

}

#endif

/* ERROR/WARNING messages:

E: Site-site interaction was not found

Internal SPPARKS error.

*/
