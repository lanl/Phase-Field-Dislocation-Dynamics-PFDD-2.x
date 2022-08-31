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

#include "mpi.h"
#include <stdio.h>
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "unistd.h"
#include "sys/stat.h"
#include "input.h"
#include "error.h"
#include "memory.h"
#include <fstream>
#include "Types.h"
#include <vector>
#include "app.h"
#include "solve.h"
#include "fft.h"
#include "random_mars.h"
#include "output.h"

#include "style_app.h"
#include "style_fft.h"
#include "style_command.h"
#include "style_diag.h"
#include "style_solve.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace PFDD_NS;

#define DELTALINE 256
#define DELTA 4

/* ---------------------------------------------------------------------- */

Input::Input(PFDD_C *pfdd_p, int argc, char **argv) : Pointers(pfdd_p)
{
  MPI_Comm_rank(world,&me);

  line = new char[MAXLINE];
  copy = new char[MAXLINE];
  work = new char[MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

  echo_screen = 0;
  echo_log = 1;

  label_active = 0;
  labelstr = NULL;
  jump_skip = 0;

  if (me == 0) {
    nfile = maxfile = 1;
    infiles = (FILE **) memory->smalloc(sizeof(FILE *),"input:infiles");
    infiles[0] = infile;
  } else infiles = NULL;

  // process command-line args
  // check for args "-var" and "-echo"
  // caller has already checked that sufficient arguments exist

  int iarg = 0;
  while (iarg < argc) {
    // if (strcmp(argv[iarg],"-var") == 0) {
    //   variable->set(argv[iarg+1],argv[iarg+2]);
    //   iarg += 3;
    // }
    if (strcmp(argv[iarg],"-echo") == 0) {
      narg = 1;
      char **tmp = arg;        // trick echo() into using argv instead of arg
      arg = &argv[iarg+1];
      echo();
      arg = tmp;
      iarg += 2;
    } else iarg++;
  }
}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  // don't free command and arg strings
  // they just point to other allocated memory

  delete [] line;
  delete [] copy;
  delete [] work;
  if (labelstr) delete [] labelstr;
  if (arg) memory->sfree(arg);
  if (infiles) memory->sfree(infiles);

}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
  int n;

  while (1) {

    // read one line from input script
    // if line ends in continuation char '&', concatenate next line(s)
    // n = str length of line

    if (me == 0) {
      if (fgets(line,MAXLINE,infile) == NULL) n = 0;
      else n = strlen(line) + 1;
      while (n >= 3 && line[n-3] == '&') {
	if (fgets(&line[n-3],MAXLINE-n+3,infile) == NULL) n = 0;
	else n = strlen(line) + 1;
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      if (label_active) error->all(FLERR,"Label wasn't found in input script");
      if (me == 0) {
	if (infile != stdin) fclose(infile);
	nfile--;
      }
      MPI_Bcast(&nfile,1,MPI_INT,0,world);
      if (nfile == 0) break;
      if (me == 0) infile = infiles[nfile-1];
      continue;
    }

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // if n = MAXLINE, line is too long

    if (n == MAXLINE) {
      char str[MAXLINE+32];
      sprintf(str,"Input line too long: %s",line);
      error->all(FLERR,str);
    }

    // echo the command unless scanning for label

    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen,"%s",line);
      if (echo_log && logfile) fprintf(logfile,"%s",line);
    }

    // parse the line
    // if no command, skip to next line in input script

    parse();
    if (command == NULL) continue;

    // if scanning for label, skip command unless it's a label command

    if (label_active && strcmp(command,"label") != 0) continue;

    // execute the command

    if (execute_command()) {
      char str[MAXLINE];
      sprintf(str,"Unknown command: %s",line);
      error->all(FLERR,str);
    }
  }
}

/* ----------------------------------------------------------------------
   process all input from filename
   called from library interface
------------------------------------------------------------------------- */

void Input::file(const char *filename)
{
  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set infile, infiles[0]

  if (me == 0) {
    if (nfile > 1)
      error->one(FLERR,"Another input script is already being processed");
    if (infile != stdin) fclose(infile);
    infile = fopen(filename,"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",filename);
      error->one(FLERR,str);
    }
    infiles[0] = infile;
  } else infile = NULL;

  file();
}

/* ----------------------------------------------------------------------
   copy command in single to line, parse and execute it
   return command name to caller
------------------------------------------------------------------------- */

char *Input::one(const char *single)
{
  strcpy(line,single);

  // echo the command unless scanning for label

  if (me == 0 && label_active == 0) {
    if (echo_screen && screen) fprintf(screen,"%s",line);
    if (echo_log && logfile) fprintf(logfile,"%s",line);
  }

  // parse the line
  // if no command, just return NULL

  parse();
  if (command == NULL) return NULL;

  // if scanning for label, skip command unless it's a label command

  if (label_active && strcmp(command,"label") != 0) return NULL;

  // execute the command and return its name

  if (execute_command()) {
    char str[MAXLINE];
    sprintf(str,"Unknown command: %s",line);
    error->all(FLERR,str);
  }
  return command;
}

/* ----------------------------------------------------------------------
   parse copy of command line
   strip comment = all chars from # on
   replace all $ via variable substitution
   command = first word
   narg = # of args
   arg[] = individual args
   treat text between double quotes as one arg
------------------------------------------------------------------------- */

void Input::parse()
{
  // make a copy to work on

  strcpy(copy,line);

  // strip any # comment by resetting string terminator
  // do not strip # inside double quotes

  int level = 0;
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#' && level == 0) {
      *ptr = '\0';
      break;
    }
    if (*ptr == '"') {
      if (level == 0) level = 1;
      else level = 0;
    }
    ptr++;
  }

  // perform $ variable substitution (print changes)
  // except if searching for a label since earlier variable may not be defined

  if (!label_active) substitute(copy,1);

  // command = 1st arg

  command = strtok(copy," \t\n\r\f");
  if (command == NULL) return;

  // point arg[] at each subsequent arg
  // treat text between double quotes as one arg
  // insert string terminators in copy to delimit args

  narg = 0;
  while (1) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = strtok(NULL," \t\n\r\f");
    if (arg[narg] && arg[narg][0] == '\"') {
      arg[narg] = &arg[narg][1];
      if (arg[narg][strlen(arg[narg])-1] == '\"')
	arg[narg][strlen(arg[narg])-1] = '\0';
      else {
	arg[narg][strlen(arg[narg])] = ' ';
	ptr = strtok(arg[narg],"\"");
	if (ptr == NULL) error->all(FLERR,"Unbalanced quotes in input line");
      }
    }
    if (arg[narg]) narg++;
    else break;
  }
}

/* ----------------------------------------------------------------------
   substitute for $ variables in str
   print updated string if flag is set and not searching for label
------------------------------------------------------------------------- */

void Input::substitute(char *str, int flag)
{
  // use work[] as scratch space to expand str
  // do not replace $ inside double quotes as flagged by level
  // var = pts at variable name, ended by NULL
  //   if $ is followed by '{', trailing '}' becomes NULL
  //   else $x becomes x followed by NULL
  // beyond = pts at text following variable

  char *var,*value,*beyond;
  int level = 0;
  char *ptr = str;

  while (*ptr) {
    if (*ptr == '$' && level == 0) {
      if (*(ptr+1) == '{') {
	var = ptr+2;
	int i = 0;
	while (var[i] != '\0' && var[i] != '}') i++;
	if (var[i] == '\0') error->one(FLERR,"Invalid variable name");
	var[i] = '\0';
	beyond = ptr + strlen(var) + 3;
      } else {
	var = ptr;
	var[0] = var[1];
	var[1] = '\0';
	beyond = ptr + strlen(var) + 1;
      }
      //value = variable->retrieve(var);
      if (value == NULL) error->one(FLERR,"Substitution for illegal variable");

      *ptr = '\0';
      strcpy(work,str);
      if (strlen(work)+strlen(value) >= MAXLINE)
	error->one(FLERR,"Input line too long after variable substitution");
      strcat(work,value);
      if (strlen(work)+strlen(beyond) >= MAXLINE)
	error->one(FLERR,"Input line too long after variable substitution");
      strcat(work,beyond);
      strcpy(str,work);
      ptr += strlen(value);
      if (flag && me == 0 && label_active == 0) {
	if (echo_screen && screen) fprintf(screen,"%s",str);
	if (echo_log && logfile) fprintf(logfile,"%s",str);
      }
      continue;
    }
    if (*ptr == '"') {
      if (level == 0) level = 1;
      else level = 0;
    }
    ptr++;
  }
}


/* ----------------------------------------------------------------------
   find next word in str
   insert 0 at end of word
   ignore leading whitespace
   treat text between single/double quotes as one arg
     matching quote must be followed by whitespace char if not end of string
     strip quotes from returned word
   return ptr to start of word
   return next = ptr after word or NULL if word ended with 0
   return NULL if no word in string
------------------------------------------------------------------------- */

char *Input::nextword(char *str, char **next)
{
  char *start,*stop;

  start = &str[strspn(str," \t\n\v\f\r")];
  if (*start == '\0') return NULL;

  if (*start == '"' || *start == '\'') {
    stop = strchr(&start[1],*start);
    if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    if (stop[1] && !isspace(stop[1]))
      error->all(FLERR,"Input line quote not followed by whitespace");
    start++;
  } else stop = &start[strcspn(start," \t\n\v\f\r")];

  if (*stop == '\0') *next = NULL;
  else *next = stop+1;
  *stop = '\0';
  return start;
}

/* ----------------------------------------------------------------------
   rellocate a string
   if n > 0: set max >= n in increments of DELTALINE
   if n = 0: just increment max by DELTALINE
------------------------------------------------------------------------- */

void Input::reallocate(char *&str, int &max, int n)
{
  if (n) {
    while (n > max) max += DELTALINE;
  } else max += DELTALINE;

  str = (char *) memory->srealloc(str,max*sizeof(char),"input:str");
}

/* ----------------------------------------------------------------------
   process a single parsed command
   return 0 if successful, -1 if did not recognize command
------------------------------------------------------------------------- */

int Input::execute_command()
{
  int flag = 1;

  if (!strcmp(command,"clear")) clear();
  else if (!strcmp(command,"echo")) echo();
  else if (!strcmp(command,"if")) ifthenelse();
  else if (!strcmp(command,"include")) include();
  else if (!strcmp(command,"jump")) jump();
  else if (!strcmp(command,"label")) label();
  else if (!strcmp(command,"log")) log();
  else if (!strcmp(command,"next")) next_command();
  else if (!strcmp(command,"print")) print();
  //else if (!strcmp(command,"variable")) variable_command();

  else if (!strcmp(command,"app_style")) app_style();
  else if (!strcmp(command,"dimension")) dimension();
  else if (!strcmp(command,"diag_style")) diag_style();
  else if (!strcmp(command,"dump")) dump();
  else if (!strcmp(command,"dump_modify")) dump_modify();
  else if (!strcmp(command,"lattice")) lattice();
  else if (!strcmp(command,"material")) material();
  else if (!strcmp(command,"fft_style")) fft_style();
  else if (!strcmp(command,"region")) region();
  else if (!strcmp(command,"run")) run();
  else if (!strcmp(command,"seed")) seed();
  else if (!strcmp(command,"solve_style")) solve_style();
  else if (!strcmp(command,"stats")) stats();

  else flag = 0;

  // return if command was listed above

  if (flag) return 0;

  // check if command is added via style.h

  if (0) return 0;      // dummy line to enable else-if macro expansion

#define COMMAND_CLASS
#define CommandStyle(key,Class)         \
  else if (strcmp(command,#key) == 0) { \
    Class key(pfdd_p);                     \
    key.command(narg,arg);              \
    return 0;                           \
  }
#include "style_command.h"
#undef COMMAND_CLASS

  // assume command is application-specific

  if (app == NULL)
    error->all(FLERR,"App_style specific command before app_style set");
  app->input(command,narg,arg);
  return 0;

  return -1;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Input::clear()
{
  if (narg > 0) error->all(FLERR,"Illegal clear command");
  pfdd_p->destroy();
  pfdd_p->create();
}

/* ---------------------------------------------------------------------- */

void Input::echo()
{
  if (narg != 1) error->all(FLERR,"Illegal echo command");

  if (strcmp(arg[0],"none") == 0) {
    echo_screen = 0;
    echo_log = 0;
  } else if (strcmp(arg[0],"screen") == 0) {
    echo_screen = 1;
    echo_log = 0;
  } else if (strcmp(arg[0],"log") == 0) {
    echo_screen = 0;
    echo_log = 1;
  } else if (strcmp(arg[0],"both") == 0) {
    echo_screen = 1;
    echo_log = 1;
  } else error->all(FLERR,"Illegal echo command");
}

/* ---------------------------------------------------------------------- */

void Input::ifthenelse()
{
  if (narg != 5 && narg != 7) error->all(FLERR,"Illegal if command");

  int flag = 0;
  if (strcmp(arg[1],"==") == 0) {
    if (atof(arg[0]) == atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],"!=") == 0) {
    if (atof(arg[0]) != atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],"<") == 0) {
    if (atof(arg[0]) < atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],"<=") == 0) {
    if (atof(arg[0]) <= atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],">") == 0) {
    if (atof(arg[0]) > atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],">=") == 0) {
    if (atof(arg[0]) >= atof(arg[2])) flag = 1;
  } else error->all(FLERR,"Illegal if command");

  if (strcmp(arg[3],"then") != 0) error->all(FLERR,"Illegal if command");
  if (narg == 7 && strcmp(arg[5],"else") != 0)
    error->all(FLERR,"Illegal if command");

  char str[128] = "\0";
  if (flag) strcpy(str,arg[4]);
  else if (narg == 7) strcpy(str,arg[6]);
  strcat(str,"\n");

  if (strlen(str) > 1) char *tmp = one(str);
}

/* ---------------------------------------------------------------------- */

void Input::include()
{
  if (narg != 1) error->all(FLERR,"Illegal include command");

  if (me == 0) {
    if (nfile == maxfile) {
      maxfile++;
      infiles = (FILE **)
        memory->srealloc(infiles,maxfile*sizeof(FILE *),"input:infiles");
    }
    infile = fopen(arg[0],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[0]);
      error->one(FLERR,str);
    }
    infiles[nfile++] = infile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::jump()
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal jump command");

  if (jump_skip) {
    jump_skip = 0;
    return;
  }

  if (me == 0) {
    if (infile != stdin) fclose(infile);
    infile = fopen(arg[0],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[0]);
      error->one(FLERR,str);
    }
    infiles[nfile-1] = infile;
  }

  if (narg == 2) {
    label_active = 1;
    if (labelstr) delete [] labelstr;
    int n = strlen(arg[1]) + 1;
    labelstr = new char[n];
    strcpy(labelstr,arg[1]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::label()
{
  if (narg != 1) error->all(FLERR,"Illegal label command");
  if (label_active && strcmp(labelstr,arg[0]) == 0) label_active = 0;
}

/* ---------------------------------------------------------------------- */

void Input::log()
{
  if (narg != 1) error->all(FLERR,"Illegal log command");

  if (me == 0) {
    if (logfile) fclose(logfile);
    if (strcmp(arg[0],"none") == 0) logfile = NULL;
    else {
      logfile = fopen(arg[0],"w");
      if (logfile == NULL) {
	char str[128];
	sprintf(str,"Cannot open logfile %s",arg[0]);
	error->one(FLERR,str);
      }
    }
    //if (universe->nworlds == 1) universe->ulogfile = logfile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::next_command()
{
  //if (variable->next(narg,arg)) jump_skip = 1;
}

/* ---------------------------------------------------------------------- */

void Input::print()
{
  if (narg != 1) error->all(FLERR,"Illegal print command");

  // substitute for $ variables (no printing)

  substitute(arg[0],0);

  if (me == 0) {
    if (screen) fprintf(screen,"%s\n",arg[0]);
    if (logfile) fprintf(logfile,"%s\n",arg[0]);
  }
}

void Input::app_style()
{
  if (fft != NULL)
    error->all(FLERR,"App_style command after simulation box is defined");
  if (narg < 1) error->all(FLERR,"Illegal app command");
  delete app;
  delete solve;
  solve = NULL;

  if (strcmp(arg[0],"none") == 0) error->all(FLERR,"Illegal app_style command");

#define APP_CLASS
#define AppStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) app = new Class(pfdd_p,narg,arg);
#include "style_app.h"
#undef APP_CLASS

  else error->all(FLERR,"Illegal app_style command");
}



/* ---------------------------------------------------------------------- */

void Input::fft_style()
{
  if (app == NULL) error->all(FLERR,"FFT_style command before app_style set");
  if (narg < 1) error->all(FLERR,"Illegal FFT_style command");
  delete fft;

  if (strcmp(arg[0],"none") == 0) fft = NULL;

#define FFT_CLASS
#define FftStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) fft = new Class(pfdd_p,narg,arg);
#include "style_fft.h"
#undef FFT_CLASS

  else error->all(FLERR,"Illegal fft_style command");
}

/* ---------------------------------------------------------------------- */

void Input::solve_style()
{
  if (app == NULL) error->all(FLERR,"Diag_style command before app_style set");

  if (narg < 1) error->all(FLERR,"Illegal solve_style command");

  if (strcmp(arg[0],"none") == 0) error->all(FLERR,"Illegal solve_style command");

#define SOLVE_CLASS
#define SolveStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) solve = new Class(pfdd_p,narg,arg);
#include "style_solve.h"
#undef SOLVE_CLASS

  else error->all(FLERR,"Illegal solve_style command");
}

void Input::diag_style()
{
  if (app == NULL) error->all(FLERR,"Diag_style command before app_style set");

  if (narg < 1) error->all(FLERR,"Illegal diag_style command");

  if (strcmp(arg[0],"none") == 0) error->all(FLERR,"Illegal diag_style command");

#define DIAG_CLASS
#define DiagStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) { \
    Diag *diagtmp = new Class(pfdd_p,narg,arg); \
    output->add_diag(diagtmp); \
  }
#include "style_diag.h"
#undef DIAG_CLASS

  else error->all(FLERR,"Illegal diag_style command");
}

/* ---------------------------------------------------------------------- */

void Input::dump()
{
  if (app == NULL) error->all(FLERR,"Dump command before app_style set");

  output->add_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump_modify()
{
  if (app == NULL) error->all(FLERR,"Dump_modify command before app_style set");

  output->dump_modify(narg,arg);
}

void Input::run()
{
  if (app == NULL) error->all(FLERR,"Run command before app_style set");

  app->run(narg,arg);
}

void Input::lattice()
{
  if (fft == NULL) error->all(FLERR,"Lattice command before fft_style set");

  fft->set_lattice(narg,arg);
}

void Input::dimension()
{
  if (fft == NULL) error->all(FLERR,"dimension command before fft_style set");

  app->set_dimension(narg,arg);
}
void Input::material()
{
  if (fft == NULL) error->all(FLERR,"Material command before fft_style set");

  fft->set_material(narg,arg);
}
void Input::boundary()
{
  if (fft->box_exist)
    error->all(FLERR,"Boundary command after simulation box is defined");
  fft->set_boundary(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::region()
{
  if (app == NULL) error->all(FLERR,"Region command before app_style set");

  fft->add_region(narg,arg);
}

void Input::seed()
{
  if (narg != 1) error->all(FLERR,"Illegal seed command");

  int seed = atoi(arg[0]);
  if (seed <= 0) error->all(FLERR,"Illegal seed command");

  ranmaster->init(seed);
}

/* ---------------------------------------------------------------------- */

void Input::stats()
{
  if (app == NULL) error->all(FLERR,"Stats command before app_style set");

  output->set_stats(narg,arg);
}
