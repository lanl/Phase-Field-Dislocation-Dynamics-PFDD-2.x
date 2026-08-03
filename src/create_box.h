/* ----------------------------------------------------------------------
PFDD -- Phase Field Dislocation Dynamics

Copyright 2021 Los Alamos National Laboratory. See LICENSE.md
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
