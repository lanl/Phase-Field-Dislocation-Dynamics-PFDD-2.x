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

#include "stdio.h"
#include "pointers.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "error.h"
#include "Types.h"

using namespace PFDD_NS;

/* ---------------------------------------------------------------------- */

Memory::Memory(PFDD_C *pfdd_p) : Pointers(pfdd_p) {}

/* ----------------------------------------------------------------------
   safe malloc
------------------------------------------------------------------------- */

void *Memory::smalloc(bigint nbytes, const char *name)
{
  if (nbytes == 0) return NULL;

  void *ptr = malloc(nbytes);

  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to allocate %d bytes for array %s",
            nbytes,name);
    error->one(FLERR,str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, bigint nbytes, const char *name)
{
  if (nbytes == 0) {
    destroy(ptr);
    return NULL;
  }

  ptr = realloc(ptr,nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to reallocate %d bytes for array %s",
            nbytes,name);
    error->one(FLERR,str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe free
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   erroneous usage of templated create/grow functions
------------------------------------------------------------------------- */

void Memory::fail(const char *name)
{
  char str[128];
  sprintf(str,"Cannot create/grow a vector/array of pointers for %s",name);
  error->one(FLERR,str);
}
