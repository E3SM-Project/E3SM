/*
** f_wrappers.c
**
** Author: Jim Rosinski
** 
** Fortran wrappers for timing library routines
*/

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <string.h>
#include <stdlib.h>
#include "private.h" /* MAX_CHARS, bool */
#include "gptl.h"    /* function prototypes */

#if ( defined FORTRANUNDERSCORE )

#define gptlinitialize gptlinitialize_
#define gptlfinalize gptlfinalize_
#define gptlpr gptlpr_
#define gptlpr_file gptlpr_file_
#define gptlpr_summary gptlpr_summary_
#define gptlpr_summary_file gptlpr_summary_file_
#define gptlbarrier gptlbarrier_
#define gptlreset gptlreset_
#define gptlstamp gptlstamp_
#define gptlstart gptlstart_
#define gptlinit_handle gptlinit_handle_
#define gptlstart_handle gptlstart_handle_
#define gptlstop gptlstop_
#define gptlstop_handle gptlstop_handle_
#define gptlsetoption gptlsetoption_
#define gptlenable gptlenable_
#define gptldisable gptldisable_
#define gptlsetutr gptlsetutr_
#define gptlquery gptlquery_
#define gptlquerycounters gptlquerycounters_
#define gptlget_wallclock gptlget_wallclock_
#define gptlget_eventvalue gptlget_eventvalue_
#define gptlget_nregions gptlget_nregions_
#define gptlget_regionname gptlget_regionname_
#define gptlget_memusage gptlget_memusage_
#define gptlprint_memusage gptlprint_memusage_
#define gptlprint_rusage gptlprint_rusage_
#define gptlnum_errors gptlnum_errors_
#define gptlnum_warn gptlnum_warn_
#define gptlget_count gptlget_count_
#define gptl_papilibraryinit gptl_papilibraryinit_
#define gptlevent_name_to_code gptlevent_name_to_code_
#define gptlevent_code_to_name gptlevent_code_to_name_

#elif ( defined FORTRANDOUBLEUNDERSCORE )

#define gptlinitialize gptlinitialize_
#define gptlfinalize gptlfinalize_
#define gptlpr gptlpr_
#define gptlpr_file gptlpr_file__
#define gptlpr_summary gptlpr_summary__
#define gptlpr_summary_file gptlpr_summary_file__
#define gptlbarrier gptlbarrier_
#define gptlreset gptlreset_
#define gptlstamp gptlstamp_
#define gptlstart gptlstart_
#define gptlinit_handle gptlinit_handle__
#define gptlstart_handle gptlstart_handle__
#define gptlstop gptlstop_
#define gptlstop_handle gptlstop_handle__
#define gptlsetoption gptlsetoption_
#define gptlenable gptlenable_
#define gptldisable gptldisable_
#define gptlsetutr gptlsetutr_
#define gptlquery gptlquery_
#define gptlquerycounters gptlquerycounters_
#define gptlget_wallclock gptlget_wallclock__
#define gptlget_eventvalue gptlget_eventvalue__
#define gptlget_nregions gptlget_nregions__
#define gptlget_regionname gptlget_regionname__
#define gptlget_memusage gptlget_memusage__
#define gptlprint_memusage gptlprint_memusage__
#define gptlprint_rusage gptlprint_rusage__
#define gptlnum_errors gptlnum_errors__
#define gptlnum_warn gptlnum_warn__
#define gptlget_count gptlget_count__
#define gptl_papilibraryinit gptl_papilibraryinit__
#define gptlevent_name_to_code gptlevent_name_to_code__
#define gptlevent_code_to_name gptlevent_code_to_name__

#endif

/* Local function prototypes */
int gptlinitialize (void);
int gptlfinalize (void);
int gptlpr (int *procid);
int gptlpr_file (char *file, int nc1);
#ifdef HAVE_MPI
int gptlpr_summary (int *fcomm);
int gptlpr_summary_file (int *fcomm, char *name, int nc1);
int gptlbarrier (int *fcomm, char *name, int nc1);
#else
int gptlpr_summary (void);
int gptlpr_summary_file (char *name, int nc1);
int gptlbarrier (void);
#endif
int gptlreset (void);
int gptlstamp (double *wall, double *usr, double *sys);
int gptlstart (char *name, int nc1);
int gptlinit_handle (char *name, int *, int nc1);
int gptlstart_handle (char *name, int *, int nc1);
int gptlstop (char *name, int nc1);
int gptlstop_handle (char *name, int *, int nc1);
int gptlsetoption (int *option, int *val);
int gptlenable (void);
int gptldisable (void);
int gptlsetutr (int *option);
int gptlquery (const char *name, int *t, int *count, int *onflg, double *wallclock, 
	       double *usr, double *sys, long long *papicounters_out, int *maxcounters, 
	       int nc);
int gptlquerycounters (const char *name, int *t, long long *papicounters_out, int nc);
int gptlget_wallclock (const char *name, int *t, double *value, int nc);
int gptlget_eventvalue (const char *timername, const char *eventname, int *t, double *value, 
			int nc1, int nc2);
int gptlget_nregions (int *t, int *nregions);
int gptlget_regionname (int *t, int *region, char *name, int nc);
int gptlget_memusage (int *size, int *rss, int *share, int *text, int *datastack);
int gptlprint_memusage (const char *str, int nc);
int gptlprint_rusage (const char *str, int nc);
int gptlnum_errors (void);
int gptlnum_warn (void);
int gptlget_count (char *, int *, int *, int);
#ifdef HAVE_PAPI
int gptl_papilibraryinit (void);
int gptlevent_name_to_code (const char *str, int *code, int nc);
int gptlevent_code_to_name (int *code, char *str, int nc);
#endif

/* Fortran wrapper functions start here */
int gptlinitialize (void)
{
  return GPTLinitialize ();
}

int gptlfinalize (void)
{
  return GPTLfinalize ();
}

int gptlpr (int *procid)
{
  return GPTLpr (*procid);
}

int gptlpr_file (char *file, int nc1)
{
  char *locfile;
  int ret;

  if ( ! (locfile = (char *) malloc (nc1+1)))
    return GPTLerror ("gptlpr_file: malloc error\n");

  snprintf (locfile, nc1+1, "%s", file);

  ret = GPTLpr_file (locfile);
  free (locfile);
  return ret;
}

#ifdef HAVE_MPI

int gptlpr_summary (int *fcomm)
{
  MPI_Comm ccomm;
#ifdef HAVE_COMM_F2C
  ccomm = MPI_Comm_f2c (*fcomm);
#else
  /* Punt and try just casting the Fortran communicator */
  ccomm = (MPI_Comm) *fcomm;
#endif
  return GPTLpr_summary (ccomm);
}

int gptlpr_summary_file (int *fcomm, char *outfile, int nc1)
{
  MPI_Comm ccomm;
  char *locfile;
  int ret;

  if ( ! (locfile = (char *) malloc (nc1+1)))
    return GPTLerror ("gptlpr_summary_file: malloc error\n");

  snprintf (locfile, nc1+1, "%s", outfile);

#ifdef HAVE_COMM_F2C
  ccomm = MPI_Comm_f2c (*fcomm);
#else
  /* Punt and try just casting the Fortran communicator */
  ccomm = (MPI_Comm) *fcomm;
#endif
  ret = GPTLpr_summary_file (ccomm, locfile);
  free (locfile);
  return ret;
}

int gptlbarrier (int *fcomm, char *name, int nc1)
{
  MPI_Comm ccomm;
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
#ifdef HAVE_COMM_F2C
  ccomm = MPI_Comm_f2c (*fcomm);
#else
  /* Punt and try just casting the Fortran communicator */
  ccomm = (MPI_Comm) *fcomm;
#endif
  return GPTLbarrier (ccomm, cname);
}

#else

int gptlpr_summary (void)
{
  return GPTLpr_summary ();
}

int gptlpr_summary_file (char *outfile, int nc1)
{
  char *locfile;
  int ret;

  if ( ! (locfile = (char *) malloc (nc1+1)))
    return GPTLerror ("gptlpr_summary_file: malloc error\n");

  snprintf (locfile, nc1+1, "%s", outfile);
  ret = GPTLpr_summary_file (locfile);
  free (locfile);
  return ret;
}

int gptlbarrier (void)
{
  return GPTLerror ("gptlbarrier: Need to build GPTL with #define HAVE_MPI to enable this routine\n");
}

#endif


int gptlreset (void)
{
  return GPTLreset ();
}

int gptlstamp (double *wall, double *usr, double *sys)
{
  return GPTLstamp (wall, usr, sys);
}

int gptlstart (char *name, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLstart (cname);
}

int gptlinit_handle (char *name, int *handle, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLinit_handle (cname, handle);
}

int gptlstart_handle (char *name, int *handle, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLstart_handle (cname, handle);
}

int gptlstop (char *name, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLstop (cname);
}

int gptlstop_handle (char *name, int *handle, int nc1)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLstop_handle (cname, handle);
}

int gptlsetoption (int *option, int *val)
{
  return GPTLsetoption (*option, *val);
}

int gptlenable (void)
{
  return GPTLenable ();
}

int gptldisable (void)
{
  return GPTLdisable ();
}

int gptlsetutr (int *option)
{
  return GPTLsetutr (*option);
}

int gptlquery (const char *name, int *t, int *count, int *onflg, double *wallclock, 
	       double *usr, double *sys, long long *papicounters_out, int *maxcounters, 
	       int nc)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLquery (cname, *t, count, onflg, wallclock, usr, sys, papicounters_out, *maxcounters);
}

int gptlquerycounters (const char *name, int *t, long long *papicounters_out, int nc)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';
  return GPTLquerycounters (cname, *t, papicounters_out);
}

int gptlget_wallclock (const char *name, int *t, double *value, int nc)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';

  return GPTLget_wallclock (cname, *t, value);
}

int gptlget_eventvalue (const char *timername, const char *eventname, int *t, double *value, 
			int nc1, int nc2)
{
  char ctimername[MAX_CHARS+1];
  char ceventname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc1, MAX_CHARS);
  strncpy (ctimername, timername, numchars);
  ctimername[numchars] = '\0';

  numchars = MIN (nc2, MAX_CHARS);
  strncpy (ceventname, eventname, numchars);
  ceventname[numchars] = '\0';

  return GPTLget_eventvalue (ctimername, ceventname, *t, value);
}

int gptlget_nregions (int *t, int *nregions)
{
  return GPTLget_nregions (*t, nregions);
}

int gptlget_regionname (int *t, int *region, char *name, int nc)
{
  int n;
  int ret;

  ret = GPTLget_regionname (*t, *region, name, nc);
  /* Turn nulls into spaces for fortran */
  for (n = 0; n < nc; ++n)
    if (name[n] == '\0')
      name[n] = ' ';
  return ret;
}

int gptlget_memusage (int *size, int *rss, int *share, int *text, int *datastack)
{
  return GPTLget_memusage (size, rss, share, text, datastack);
}

int gptlprint_memusage (const char *str, int nc)
{
  char cname[128+1];
  int numchars = MIN (nc, 128);

  strncpy (cname, str, numchars);
  cname[numchars] = '\0';
  return GPTLprint_memusage (cname);
}

int gptlprint_rusage (const char *str, int nc)
{
  char cname[128+1];
  int numchars = MIN (nc, 128);

  strncpy (cname, str, numchars);
  cname[numchars] = '\0';
  return GPTLprint_rusage (cname);
}

int gptlnum_errors (void)
{
  return GPTLnum_errors ();
}

int gptlnum_warn (void)
{
  return GPTLnum_warn ();
}

int gptlget_count (char *name, int *t, int *count, int nc)
{
  char cname[MAX_CHARS+1];
  int numchars;

  numchars = MIN (nc, MAX_CHARS);
  strncpy (cname, name, numchars);
  cname[numchars] = '\0';

  return GPTLget_count (cname, *t, count);
}

#ifdef HAVE_PAPI
#include <papi.h>

int gptl_papilibraryinit (void)
{
  return GPTL_PAPIlibraryinit ();;
}

int gptlevent_name_to_code (const char *str, int *code, int nc)
{
  char cname[PAPI_MAX_STR_LEN+1];
  int numchars = MIN (nc, PAPI_MAX_STR_LEN);

  strncpy (cname, str, numchars);
  cname[numchars] = '\0';

  /* "code" is an int* and is an output variable */
  return GPTLevent_name_to_code (cname, code);
}

int gptlevent_code_to_name (int *code, char *str, int nc)
{
  int i;

  if (nc < PAPI_MAX_STR_LEN)
    return GPTLerror ("gptl_event_code_to_name: output name must hold at least %d characters\n",
		      PAPI_MAX_STR_LEN);

  if (GPTLevent_code_to_name (*code, str) == 0) {
    for (i = strlen(str); i < nc; ++i)
      str[i] = ' ';
  } else {
    return GPTLerror ("");
  }
  return 0;
}
#endif
