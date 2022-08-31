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

#ifndef PFDD_INPUT_H
#define PFDD_INPUT_H

#include "stdio.h"
#include "pointers.h"
#include "mpi.h"

#include <vector>
#include "error.h"

namespace PFDD_NS {

  class Input : protected Pointers {
  public:
    int narg;                    // # of command args
    char **arg;                  // parsed args for command
    FILE **infiles;

    Input(class PFDD_C *, int, char **);
    ~Input();
    void file();                   // process all input
    void file(const char *);       // process an input script
    char *one(const char *);       // process a single command
    void substitute(char *, int);
    // substitute for variables in a string

  private:
    int me;                      // proc ID
    char *command;               // ptr to current command
    int maxarg;                  // max # of args in arg
    char *line,*copy,*work;      // input line & copy and work string
    int echo_screen;             // 0 = no, 1 = yes
    int echo_log;                // 0 = no, 1 = yes
    int nfile,maxfile;           // current # and max # of open input files
    int label_active;            // 0 = no label, 1 = looking for label
    char *labelstr;              // label string being looked for
    int jump_skip;               // 1 if skipping next jump, 0 otherwise
    int gamma_flag;              // 1 if GammaDotBase already read
    int rate_flag;               // 1 if RateSensitivity already read
    int ready_flag;              // 1 if dimension can be called

    void parse();                          // parse an input text line
    char *nextword(char *, char **);       // find next word in string with quotes
    void reallocate(char *&, int &, int);  // reallocate a char string
    int execute_command();                 // execute a single command

    void clear();                // individual commands
    void echo();
    void ifthenelse();
    void include();
    void jump();
    void label();
    void log();
    void next_command();
    void print();
    //void variable_command();

    void app_style();
    void lattice();
    void diag_style();
    void dimension();
    void dump();
    void dump_modify();
    void boundary();
    void material();
    void region();
    void fft_style();
    void solve_style();
    void run();
    void seed();
    void stats();

  };

}

#endif

/* ERROR/WARNING messages:

E: Label wasn't found in input script

Self-explanatory.

E: Unknown command: %s

The command is not known to LAMMPS.  Check the input script.

E: Another input script is already being processed

Cannot attempt to open a 2nd input script, when the original file is
still being processed.

E: Cannot open input script %s

Self-explanatory.

E: Unbalanced quotes in input line

No matching end double quote was found following a leading double
quote.

E: Input line quote not followed by whitespace

An end quote must be followed by whitespace.

E: Invalid variable name

Variable name used in an input script line is invalid.

E: Invalid immediate variable

Syntax of immediate value is incorrect.

E: Substitution for illegal variable

Input script line contained a variable that could not be substituted
for.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open logfile %s

The LAMMPS log file specified in the input script cannot be opened.
Check that the path and name are correct.

E: Angle_coeff command before simulation box is defined

The angle_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Angle_coeff command before angle_style is defined

Coefficients cannot be set in the data file or via the angle_coeff
command until an angle_style has been assigned.

E: Angle_coeff command when no angles allowed

The chosen atom style does not allow for angles to be defined.

E: Angle_style command when no angles allowed

The chosen atom style does not allow for angles to be defined.

E: Atom_style command after simulation box is defined

The atom_style command cannot be used after a read_data,
read_restart, or create_box command.

E: Bond_coeff command before simulation box is defined

The bond_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Bond_coeff command before bond_style is defined

Coefficients cannot be set in the data file or via the bond_coeff
command until an bond_style has been assigned.

E: Bond_coeff command when no bonds allowed

The chosen atom style does not allow for bonds to be defined.

E: Bond_style command when no bonds allowed

The chosen atom style does not allow for bonds to be defined.

E: Boundary command after simulation box is defined

The boundary command cannot be used after a read_data, read_restart,
or create_box command.

E: Box command after simulation box is defined

The box command cannot be used after a read_data, read_restart, or
create_box command.

E: Dihedral_coeff command before simulation box is defined

The dihedral_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Dihedral_coeff command before dihedral_style is defined

Coefficients cannot be set in the data file or via the dihedral_coeff
command until an dihedral_style has been assigned.

E: Dihedral_coeff command when no dihedrals allowed

The chosen atom style does not allow for dihedrals to be defined.

E: Dihedral_style command when no dihedrals allowed

The chosen atom style does not allow for dihedrals to be defined.

E: Dimension command after simulation box is defined

The dimension command cannot be used after a read_data,
read_restart, or create_box command.

E: Improper_coeff command before simulation box is defined

The improper_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Improper_coeff command before improper_style is defined

Coefficients cannot be set in the data file or via the improper_coeff
command until an improper_style has been assigned.

E: Improper_coeff command when no impropers allowed

The chosen atom style does not allow for impropers to be defined.

E: Improper_style command when no impropers allowed

The chosen atom style does not allow for impropers to be defined.

E: KSpace style has not yet been set

Cannot use kspace_modify command until a kspace style is set.

E: Mass command before simulation box is defined

The mass command cannot be used before a read_data, read_restart, or
create_box command.

E: Min_style command before simulation box is defined

The min_style command cannot be used before a read_data, read_restart,
or create_box command.

E: Newton bond change after simulation box is defined

The newton command cannot be used to change the newton bond value
after a read_data, read_restart, or create_box command.

E: Package command after simulation box is defined

The package command cannot be used afer a read_data, read_restart, or
create_box command.

E: Package cuda command without USER-CUDA installed

The USER-CUDA package must be installed via "make yes-user-cuda"
before LAMMPS is built.

E: Pair_coeff command before simulation box is defined

The pair_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Pair_coeff command before pair_style is defined

Self-explanatory.

E: Pair_modify command before pair_style is defined

Self-explanatory.

E: Pair_write command before pair_style is defined

Self-explanatory.

E: Processors command after simulation box is defined

The processors command cannot be used after a read_data, read_restart,
or create_box command.

E: Run_style command before simulation box is defined

The run_style command cannot be used before a read_data,
read_restart, or create_box command.

E: Units command after simulation box is defined

The units command cannot be used after a read_data, read_restart, or
create_box command.

*/
