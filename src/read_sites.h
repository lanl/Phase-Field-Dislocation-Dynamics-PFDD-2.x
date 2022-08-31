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
CommandStyle(read_sites,ReadSites)

#else

#ifndef PFDD_READ_SITES_H
#define PFDD_READ_SITES_H

#include "pointers.h"

namespace PFDD_NS {

class ReadSites : protected Pointers {
 public:
  ReadSites(class PFDD_C *);
  ~ReadSites();
  void command(int, char **);

 private:
  int me;
  char *line,*keyword,*buffer;
  FILE *fp;
  int narg,maxarg,compressed;
  char **arg;

  int latticeflag;

  int maxneigh;
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;

  void open(char *);
  void header();
  void parse_keyword(int);
  void parse_coeffs(int, char *);

  void sites();
  void neighbors();
  void values();

  int count_words(char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Read_sites command before app_style set

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot run 2d simulation with nonperiodic Z dimension

UNDOCUMENTED

E: Cannot run 1d simulation with nonperiodic Y or Z dimension

UNDOCUMENTED

E: Cannot read Sites after sites already exist

UNDOCUMENTED

E: Cannot read Neighbors after sites already exist

UNDOCUMENTED

E: Can only read Neighbors for on-lattice applications

UNDOCUMENTED

E: Cannot read Neighbors unless max neighbors is set

UNDOCUMENTED

E: Must read Sites before Neighbors

Self-explanatory.

E: Cannot read Values before sites exist or are read

UNDOCUMENTED

E: Unknown identifier in data file: %s

Self-explanatory.

E: Site file has no Sites, Neighbors, or Values

UNDOCUMENTED

E: No Sites defined in site file

UNDOCUMENTED

E: No Neighbors defined in site file

UNDOCUMENTED

E: Unexpected end of data file

Self-explanatory.

E: Data file dimension does not match existing box

UNDOCUMENTED

E: Data file number of sites does not match existing sites

UNDOCUMENTED

E: Off-lattice application data file cannot have maxneigh setting

UNDOCUMENTED

E: Data file maxneigh setting does not match existing sites

UNDOCUMENTED

E: Data file simluation box different that current box

UNDOCUMENTED

E: System in site file is too big

UNDOCUMENTED

E: Incorrect site format in data file

Self-explanatory.

E: Did not assign all sites correctly

One or more sites in the read_sites file were not assigned to 
a processor correctly.

E: Invalid site ID in Sites section of data file

Self-explanatory.

E: Too many neighbors per site

Internal SPPARKS error.

E: Incorrect value format in data file

Self-explanatory.

E: Cannot open gzipped file

Self-explantory.

E: Cannot open file %s

Self-explanatory.

*/
