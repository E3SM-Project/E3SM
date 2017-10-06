/*
** Print to stdout a dependency list for input file specified on the command
** line.  A dependency is anything that is referenced by a "#include"' or
** f90 "use" statement.  In addition to these dependencies, write a dependency
** rule of "file.d" for each "file.F" or "file.c".  This is to accomodate the
** default "make" procedure for CCM.
**
** The name of the module being "use"d is assumed to be case sensitive even
** though the Fortran language is not.  In addition, Fortran source files are
** assumed to end in .F.  For example, the statement "use Xxx" will translate
** into a dependency of Xxx.o, and the file searched for will be Xxx.F.
**
** Only files which exist in at least one directory named in the current
** directory or one or more "-I" command line arguments will be considered.
**
** An ANSI C compiler is required to build this code.
*/

#include <stdio.h>         /* printf, puts */
#include <stdlib.h>        /* malloc, getopt */
#include <string.h>        /* strcpy */
#include <unistd.h>        /* access */
#include <ctype.h>         /* isspace, isalnum, tolower */

#define MAXLEN 256
#define TRUE 1
#define FALSE 0

/*
** Linked list struct used for directories to search, and filenames already
** found.
*/

struct node {
  char *name;
  struct node *next;
};

/*
** lists of dependencies already found: prevents duplicates.
*/

static struct node *list = NULL;     /* For #include */
static struct node *uselist = NULL;  /* For use */
static struct node *suffix_list;      /* List of Fortran suffixes to look for */

/*
** Function prototypes
*/

static void check (char *, struct node *, char *, int);
static int already_found (char *, struct node *);

main (int argc, char **argv)
{
  int lastdot;             /* points to the last . in fname */
  int c;                   /* return from getopt */
  int recursive = FALSE;   /* flag asks for recursive check:
			   ** i.e. check the thing being #included for #includes */
  FILE *fpFname;

  char line[MAXLEN];         /* line read from input file */
  char doto[MAXLEN];         /* name of .o file (from input file) */
  char dotd[MAXLEN];         /* name of .o file (from input file) */
  char fullpath[MAXLEN];     /* full pathname to potential dependency for .F files */
  char depnam[MAXLEN];       /* dependency name (from #include or use) */
  char srcfile[MAXLEN];      /* source file .F   name from "use" */
  char *lptr;                /* points into line */
  char *fname;               /* input file name from command line */
  char *fptr;                /* pointer to copy into depnam */
  char *relpath;             /* input file name or path to it */

  struct node *dirlist;    /* list of directories to search */
  struct node *dirptr;     /* loop through dirlist */
  struct node *newnode;    /* malloc'd node */
  struct node *last;       /* last entry in #include list of found dependencies */
  struct node *uselast;    /* last entry in "use" list of found dependencies */
  struct node *sptr;       /* pointer into suffix_list */

  /*
  ** Always put "." first in Filepath since gnumake will put "." first
  ** regardless of whether it is specified in VPATH
  */

  dirlist = dirptr = (struct node *) malloc (sizeof (struct node));
  dirptr->name = (char *) malloc (2);
  strcpy (dirptr->name, ".");
  dirptr->next = NULL;

  /*
  ** Always look for .F and .F90 files.  List can be augmented via "-s" cmd line arg(s).
  */

  suffix_list = (struct node *) malloc (sizeof (struct node));
  suffix_list-> name = (char *) malloc (3);
  strcpy (suffix_list->name, ".F");

  suffix_list->next = (struct node *) malloc (sizeof (struct node));
  sptr = suffix_list->next;
  sptr->name = (char *) malloc (5);
  strcpy (sptr->name, ".F90");
  sptr->next = NULL;

  while ((c = getopt (argc, argv, "I:rs:f")) != -1) {

    switch(c) {

    case 'f':   /* this arg is for backward compatibility */
      break;
    case 'I':
      dirptr->next = (struct node *) malloc (sizeof (struct node));
      dirptr = dirptr->next;
      dirptr->name = (char *) malloc (strlen (optarg) + 1);
      strcpy (dirptr->name, optarg);
      dirptr->next = NULL;
      break;
    case 's':
      sptr->next = (struct node *) malloc (sizeof (struct node));
      sptr = sptr->next;
      sptr->name = (char *) malloc (strlen (optarg) + 2);
      strcpy (sptr->name, ".");
      strcat (sptr->name, optarg);
      sptr->next = NULL;
      break;
    case 'r':
      recursive = TRUE;
      break;
    case '?':   /* Unknown option */
      fprintf (stderr, "%s: Unknown option encountered\n", argv[0]);
    }
  }

  if (argc == optind+1) {
    relpath = argv[optind];

  } else {

    fprintf (stderr, "Usage: %s [-Idir] [-r] [-s suffix] file\n", argv[0]);
    exit (-1);
  }

  /*
  ** Retain only the filename of the input file for which dependencies are
  ** being generated.
  */

  fname = relpath + strlen (relpath) - 1;
  while (*fname != '/' && fname > relpath) fname--;
  if (*fname == '/') fname++;

  /*
  ** Define the .o file by changing tail to ".o"
  */

  strcpy (doto, fname);
  for (lastdot = strlen (fname) - 1; doto[lastdot] != '.' && lastdot > 0;
       lastdot--);

  if (lastdot == 0) {
    fprintf (stderr, "Input file %s needs a head\n", fname);
    exit (1);
  }

  doto[lastdot] = '\0';
  strcpy (dotd, doto);
  strcat (doto, ".o ");
  strcat (dotd, ".d ");

  /*
  ** write the blah.o blah.d: blah.F (or .c or whatever) dependency to stdout
  */

  fputs (doto   , stdout);
  fputs (dotd   , stdout);
  fputs (": "   , stdout);
  fputs (fname  , stdout);
  fputs ("\n"   , stdout);

  if ((fpFname = fopen (relpath, "r")) == NULL) {
    fprintf (stderr, "Can't open file %s\n", relpath);
    exit (1);
  }

  while (fgets (line, MAXLEN, fpFname) != NULL) {

    /*
    ** Check for dependencies of the cpp "include" variety.  Allow for lines
    ** of the form "# include"
    */

    if (line[0] == '#') {
      for (lptr = line+1; isspace (*lptr); lptr++);
      if (strncmp (lptr, "include ", 8) == 0) {
	for (lptr += 8; *lptr != '<' && *lptr != '"' && *lptr != '\0'; lptr++);

	if (*lptr == '\0')
	  break;              /* Bad input line: ignore */

	/*
	** Fill in depnam with the dependency (i.e. the thing being
	** #included.  Syntax check is not perfect.
	*/

	for (fptr = depnam; *++lptr != '>' && *lptr != '"' && *lptr != '\0';
	     fptr++)
	  *fptr = *lptr;

	if (*lptr == '\0')
	  break;              /* Bad input line: ignore */

	*fptr = '\0';

	if ( ! already_found (depnam, list)) { /* Skip any duplicates */

	  /*
	  ** Include only dependencies which are specified by -Ixxx on the
	  ** command line.  These directories are defined by the linked list
	  ** pointed to by dirlist.
	  */

	  for (dirptr = dirlist; dirptr != NULL; dirptr = dirptr->next) {
	    strcpy (fullpath, dirptr->name);
	    strcat (fullpath, "/");
	    strcat (fullpath, depnam);

	    /*
	    ** If the file exists and is readable, add an entry to the "found"
            ** list, then write a dependency rule to stdout.
	    */

	    if (access (fullpath, R_OK) == 0) {
	      newnode = malloc (sizeof (struct node));
	      newnode->name = malloc (strlen (depnam) + 1);
	      strcpy (newnode->name, depnam);
	      newnode->next = NULL;

	      if (list == NULL)
		list = newnode;
	      else
		last->next = newnode;

	      last = newnode;
	      fputs (doto , stdout);
	      fputs (": " , stdout);
	      fputs (depnam, stdout);
	      fputs ("\n", stdout);

	      /*
	      ** Check for nested #include's if flag was set
	      */

	      if (recursive) check (fullpath, dirlist, doto, 0);

	      break;  /* Dependency found: process next line */
	    }
	  }
	}
      }

    } else {

      /*
      ** Check for dependencies of the f90 "use" variety.  To strictly adhere
      ** to fortran std, should allow for spaces between chars of "use".
      */

      for (lptr = line; isspace (*lptr); lptr++);
      if (tolower ((int) lptr[0]) == 'u' &&
	  tolower ((int) lptr[1]) == 's' &&
	  tolower ((int) lptr[2]) == 'e') {

	for (lptr += 3; isspace (*lptr); lptr++);

	/*
	** Fill in depnam with the dependency (i.e. the thing being "use"d.
	** Strictly speaking, should disallow numeric starting character.
	*/

	for (fptr = depnam; isalnum (*lptr) || *lptr == '_'; (fptr++, lptr++))
	  *fptr = *lptr;
	*fptr = '\0';

	/*
	** srcfile is the source file name from which the dependency is
	** generated.  Note case sensitivity of depnam.
	*/

	if ( ! already_found (depnam, uselist)) {  /* Skip any duplicates */

	  /*
	  ** Loop through suffix list
	  */

	  for (sptr = suffix_list; sptr != NULL; sptr = sptr->next) {

	    strcpy (srcfile, depnam);
	    strcat (srcfile, sptr->name);

	    /*
	    ** Include only dependencies which are specified by -Ixxx on the
	    ** command line.  These directories are defined by the linked list
	    ** pointed to by dirlist.
	    */

	    for (dirptr = dirlist; dirptr != NULL; dirptr = dirptr->next) {
	      strcpy (fullpath, dirptr->name);
	      strcat (fullpath, "/");
	      strcat (fullpath, srcfile);

	      /*
	      ** If the file exists and is readable, add an entry to the "found"
	      ** list, then write a dependency rule to stdout.
	      */

	      if (access (fullpath, R_OK) == 0) {
		newnode = malloc (sizeof (struct node));
		newnode->name = malloc (strlen (srcfile) + 1);
		strcpy (newnode->name, depnam);
		newnode->next = NULL;

		if (uselist == NULL)
		  uselist = newnode;
		else
		  uselast->next = newnode;

		uselast = newnode;

		fputs (doto  , stdout);
		fputs (": "  , stdout);
		fputs (depnam, stdout);
		fputs (".o"  , stdout);
		fputs ("\n"  , stdout);

		goto read_next_line;  /* Dependency found: process next line */

	      } /* if (access (fullpath... */

	    }   /* loop through linked list of directories from Filepath */
	  }     /* loop through linked list of suffixes */
	}       /* if ( ! already_found (srcfile... */
      }         /* if (lptr points to "use " */
    }           /* else branch of if (line[0] == '#') */
  read_next_line:
    continue;
  }             /* Looping over lines in the file */

  fclose (fpFname);
  return (0);
}

void check (char *file, struct node *dirlist, char *doto, int recurse_level)
{
  FILE *fpFile;

  char line[MAXLEN], fullpath[MAXLEN];
  char depnam[MAXLEN];
  char *lptr, *fptr;

  struct node *dirptr;

  /*
  ** Don't bother checking beyond 3 levels of recursion
  */

  if (recurse_level > 3) {
    fprintf (stderr, "More than 3 levels of recursion detected: bailing out\n");
    return;
  }

  if ((fpFile = fopen (file, "r")) == NULL) {
    fprintf (stderr, "Can't open file %s\n", file);
    exit (1);
  }

  while (fgets (line, MAXLEN, fpFile) != NULL) {

    /*
    ** Check for dependencies of the cpp "include" variety.  Allow for lines
    ** of the form "# include"
    */

    if (line[0] == '#') {
      for (lptr = line+1; isspace (*lptr); lptr++);
      if (strncmp (lptr, "include ", 8) == 0) {
	for (lptr += 8; *lptr != '<' && *lptr != '"' && *lptr != '\0'; lptr++);

	if (*lptr == '\0')
	  break;              /* Bad input line: ignore */

	/*
	** Fill in depnam with the dependency (i.e. the thing being
	** #included.  Syntax check is not perfect.
	*/

	for (fptr = depnam; *++lptr != '>' && *lptr != '"' && *lptr != '\0'; fptr++)
	  *fptr = *lptr;

	if (*lptr == '\0')
	  break;              /* Bad input line: ignore */

	*fptr = '\0';

	/*
        ** Don't include dependencies which are not in the Filepath
        */

	for (dirptr = dirlist; dirptr != NULL; dirptr = dirptr->next) {
	  strcpy (fullpath, dirptr->name);
	  strcat (fullpath, "/");
	  strcat (fullpath, depnam);

	  /*
	  ** If the file exists and is readable, add an entry to the "found"
	  ** list, then write a dependency rule to stdout.
	  */

	  if (access (fullpath, R_OK) == 0) {
	    fputs (doto , stdout);
	    fputs (": " , stdout);
	    fputs (depnam, stdout);
	    fputs ("\n", stdout);

	    /*
	    ** Check for nested #include's
	    */

	    check (fullpath, dirlist, doto, recurse_level+1);
	    break;
	  }
	}
      }
    }
  }
  fclose (fpFile);
  return;
}

int already_found (char *name, struct node *list)
{
  struct node *ptr;

  for (ptr = list; ptr != NULL; ptr = ptr->next) {
    if (strcmp (ptr->name, name) == 0) return (1);
  }
  return (0);
}
