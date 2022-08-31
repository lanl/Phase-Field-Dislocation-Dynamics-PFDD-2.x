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

#ifdef COMMAND_CLASS
CommandStyle(create_box,CreateBox)

#else

#ifndef PFDD_CREATE_BOX_H
#define PFDD_CREATE_BOX_H

#include "pointers.h"

namespace PFDD_NS {

class CreateBox : protected Pointers {
 public:
  CreateBox(class PFDD_C *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Create_box command before app_style set

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot create box with this application style

This application does not support spatial domains.

E: Cannot create box after simulation box is defined

Self-explanatory.

E: Cannot run 2d simulation with nonperiodic Z dimension

UNDOCUMENTED

E: Cannot run 1d simulation with nonperiodic Y or Z dimension

UNDOCUMENTED

E: Create_box region ID does not exist

Self-explanatory.

E: Create_box region must be of type inside

Self-explanatory.

*/
