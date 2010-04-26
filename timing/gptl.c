/*
** $Id: gptl.c,v 1.153 2010/04/13 21:41:50 rosinski Exp $
**
** Author: Jim Rosinski
**
** Main file contains most user-accessible GPTL functions
*/
#ifndef inline
#define inline
#endif
 
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <stdlib.h>        /* malloc */
#include <sys/time.h>      /* gettimeofday */
#include <sys/times.h>     /* times */
#include <unistd.h>        /* gettimeofday, syscall */
#include <stdio.h>
#include <string.h>        /* memset, strcmp (via STRMATCH) */
#include <ctype.h>         /* isdigit */
#include <sys/types.h>     /* u_int8_t, u_int16_t */
#include <assert.h>

#ifdef HAVE_PAPI
#include <papi.h>          /* PAPI_get_real_usec */
#endif

#ifdef UNICOSMP
#include <intrinsics.h>    /* rtc */
#endif

#ifdef HAVE_LIBRT
#include <time.h>
#endif

#include "private.h"
#include "gptl.h"

static Timer **timers = 0;          /* linked list of timers */
static Timer **last = 0;            /* last element in list */
static int *max_depth;              /* maximum indentation level encountered */
static int *max_name_len;           /* max length of timer name */
static int GPTLnthreads= -1;        /* num threads. Init to bad value */
static int maxthreads  = -1;        /* max threads (=GPTLnthreads for OMP). Init to bad value */
static int depthlimit  = 99999;     /* max depth for timers (99999 is effectively infinite) */
static volatile bool disabled = false;           /* Timers disabled? */
static volatile bool initialized = false;        /* GPTLinitialize has been called */
static volatile bool pr_has_been_called = false; /* GPTLpr_file has been called */
static Entry eventlist[MAX_AUX];    /* list of PAPI-based events to be counted */
static int nevents = 0;             /* number of PAPI events (init to 0) */
static bool dousepapi = false;      /* saves a function call if stays false */
static bool verbose = false;        /* output verbosity */
static bool percent = false;        /* print wallclock also as percent of 1st timers[0] */
static bool dopr_preamble = true;   /* whether to print preamble info */
static bool dopr_threadsort = true; /* whether to print sorted thread stats */
static bool dopr_multparent = true; /* whether to print multiple parent info */
static bool dopr_collision = true;  /* whether to print hash collision info */

static time_t ref_gettimeofday = -1; /* ref start point for gettimeofday */
static time_t ref_clock_gettime = -1;/* ref start point for clock_gettime */
static long long ref_papitime = -1;  /* ref start point for PAPI_get_real_usec */

typedef struct {
  const Option option;  /* wall, cpu, etc. */
  const char *str;      /* descriptive string for printing */
  bool enabled;         /* flag */
} Settings;

/* For Summary stats */

typedef struct {
  double wallmax;
  double wallmin;
#ifdef HAVE_PAPI
  double papimax[MAX_AUX];
  double papimin[MAX_AUX];
#endif
  unsigned long count;
  int wallmax_p;               /* over processes */
  int wallmax_t;               /* over threads */
  int wallmin_p;
  int wallmin_t;
#ifdef HAVE_PAPI
  int papimax_p[MAX_AUX];      /* over processes */
  int papimax_t[MAX_AUX];      /* over threads */
  int papimin_p[MAX_AUX];
  int papimin_t[MAX_AUX];
#endif
} Summarystats;

/* Options, print strings, and default enable flags */

static Settings cpustats =      {GPTLcpu,      "Usr       sys       usr+sys   ", false};
static Settings wallstats =     {GPTLwall,     "Wallclock max       min       ", true };
static Settings overheadstats = {GPTLoverhead, "UTR_Overhead  "                , true };

static Hashentry **hashtable;    /* table of entries */
static long ticks_per_sec;       /* clock ticks per second */

typedef struct {
  int val;                       /* depth in calling tree */
  int padding[31];               /* padding is to mitigate false cache sharing */
} Nofalse; 
static Timer ***callstack;       /* call stack */
static Nofalse *stackidx;        /* index into callstack: */

static Method method = GPTLmost_frequent;  /* default parent/child printing mechanism */

/* Local function prototypes */

static void printstats (const Timer *, FILE *, const int, const int, const bool, double);
static void add (Timer *, const Timer *);

#ifdef HAVE_MPI
static void get_threadstats (const char *, Summarystats *);
static void get_summarystats (Summarystats *, const Summarystats *, int);
#endif

static void print_multparentinfo (FILE *, Timer *);
static inline int get_cpustamp (long *, long *);
static int newchild (Timer *, Timer *);
static int get_max_depth (const Timer *, const int);
static int num_descendants (Timer *);
static int is_descendant (const Timer *, const Timer *);
static char *methodstr (Method);

/* These are the (possibly) supported underlying wallclock timers */

static inline double utr_nanotime (void);
static inline double utr_rtc (void);
static inline double utr_mpiwtime (void);
static inline double utr_clock_gettime (void);
static inline double utr_papitime (void);
static inline double utr_gettimeofday (void);

static int init_nanotime (void);
static int init_rtc (void);
static int init_mpiwtime (void);
static int init_clock_gettime (void);
static int init_papitime (void);
static int init_gettimeofday (void);

static double utr_getoverhead (void);
static inline Timer *getentry_instr (const Hashentry *, void *, int *);
static inline Timer *getentry (const Hashentry *, const char *, int *);
static void printself_andchildren (const Timer *, FILE *, const int, const int, const double);
static inline int update_parent_info (Timer *, Timer **, int);
static inline int update_stats (Timer *, const double, const long, const long, const int);
static inline int update_ll_hash (Timer *, const int, const int);
static inline int update_ptr (Timer *, const int);
static int construct_tree (Timer *, Method);
static int get_max_depth (const Timer *, const int);

typedef struct {
  const Funcoption option;
  double (*func)(void);
  int (*funcinit)(void);
  const char *name;
} Funcentry;

static Funcentry funclist[] = {
  {GPTLgettimeofday, utr_gettimeofday,  init_gettimeofday,  "gettimeofday"},
  {GPTLnanotime,     utr_nanotime,      init_nanotime,      "nanotime"},
  {GPTLrtc,          utr_rtc,           init_rtc,           "_rtc"},
  {GPTLmpiwtime,     utr_mpiwtime,      init_mpiwtime,      "MPI_Wtime"},
  {GPTLclockgettime, utr_clock_gettime, init_clock_gettime, "clock_gettime"},
  {GPTLpapitime,     utr_papitime,      init_papitime,      "PAPI_get_real_usec"}
};
static const int nfuncentries = sizeof (funclist) / sizeof (Funcentry);

static double (*ptr2wtimefunc)() = 0;             /* init to invalid */
static int funcidx = 0;                           /* default timer is gettimeofday*/  

#ifdef HAVE_NANOTIME
static float cpumhz = -1.;                        /* init to bad value */
static double cyc2sec = -1;                       /* init to bad value */
static unsigned inline long long nanotime (void); /* read counter (assembler) */
static float get_clockfreq (void);                /* cycles/sec */
#endif

#ifdef UNICOSMP
static double ticks2sec = -1;                     /* init to bad value */
#endif

static const int tablesize = 128*MAX_CHARS;       /* 128 is size of ASCII char set */
static char *outdir = 0;                          /* dir to write output files to */

/*
** GPTLsetoption: set option value to true or false.
**
** Input arguments:
**   option: option to be set
**   val:    value to which option should be set (nonzero=true, zero=false)
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLsetoption (const int option,  /* option */
		   const int val)     /* value */
{
  if (initialized)
    return GPTLerror ("GPTLsetoption: must be called BEFORE GPTLinitialize\n");

  if (option == GPTLabort_on_error) {
    GPTLset_abort_on_error ((bool) val);
    if (verbose)
      printf ("GPTLsetoption: boolean abort_on_error = %d\n", val);
    return 0;
  }

  switch (option) {
  case GPTLcpu:
#ifdef HAVE_TIMES
    cpustats.enabled = (bool) val; 
    if (verbose)
      printf ("GPTLsetoption: cpustats = %d\n", val);
#else
    if (val)
      return GPTLerror ("GPTLsetoption: times() not available\n");
#endif
    return 0;
  case GPTLwall:     
    wallstats.enabled = (bool) val; 
    if (verbose)
      printf ("GPTLsetoption: boolean wallstats = %d\n", val);
    return 0;
  case GPTLoverhead: 
    overheadstats.enabled = (bool) val; 
    if (verbose)
      printf ("GPTLsetoption: boolean overheadstats = %d\n", val);
    return 0;
  case GPTLdepthlimit: 
    depthlimit = val; 
    if (verbose)
      printf ("GPTLsetoption: depthlimit = %d\n", val);
    return 0;
  case GPTLverbose: 
    verbose = (bool) val; 
#ifdef HAVE_PAPI
    (void) GPTL_PAPIsetoption (GPTLverbose, val);
#endif
    if (verbose)
      printf ("GPTLsetoption: boolean verbose = %d\n", val);
    return 0;
  case GPTLpercent: 
    percent = (bool) val; 
    if (verbose)
      printf ("GPTLsetoption: boolean percent = %d\n", val);
    return 0;
  case GPTLdopr_preamble: 
    dopr_preamble = (bool) val; 
    if (verbose)
      printf ("GPTLsetoption: boolean dopr_preamble = %d\n", val);
    return 0;
  case GPTLdopr_threadsort: 
    dopr_threadsort = (bool) val; 
    if (verbose)
      printf ("GPTLsetoption: boolean dopr_threadsort = %d\n", val);
    return 0;
  case GPTLdopr_multparent: 
    dopr_multparent = (bool) val; 
    if (verbose)
      printf ("GPTLsetoption: boolean dopr_multparent = %d\n", val);
    return 0;
  case GPTLdopr_collision: 
    dopr_collision = (bool) val; 
    if (verbose)
      printf ("GPTLsetoption: boolean dopr_collision = %d\n", val);
    return 0;
  case GPTLprint_method:
    method = (Method) val; 
    if (verbose)
      printf ("GPTLsetoption: print_method = %s\n", methodstr (method));
    return 0;
  case GPTLsync_mpi:
#ifdef ENABLE_PMPI
    if (GPTLpmpi_setoption (option, val) != 0)
      fprintf (stderr, "GPTLsetoption: GPTLpmpi_setoption failure\n");
#endif
    if (verbose)
      printf ("GPTLsetoption: boolean sync_mpi = %d\n", val);
    return 0;

  /* 
  ** Allow GPTLmultiplex to fall through because it will be handled by 
  ** GPTL_PAPIsetoption()
  */

  case GPTLmultiplex:
  default:
    break;
  }

#ifdef HAVE_PAPI
  if (GPTL_PAPIsetoption (option, val) == 0) {
    if (val)
      dousepapi = true;
    return 0;
  }
#else
  /* Make GPTLnarrowprint a placebo if PAPI not enabled */

  if (option == GPTLnarrowprint)
    return 0;
#endif

  return GPTLerror ("GPTLsetoption: faiure to enable option %d\n", option);
}

/*
** GPTLsetutr: set underlying timing routine.
**
** Input arguments:
**   option: index which sets function
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLsetutr (const int option)
{
  int i;

  if (initialized)
    return GPTLerror ("GPTLsetutr: must be called BEFORE GPTLinitialize\n");

  for (i = 0; i < nfuncentries; i++) {
    if (option == (int) funclist[i].option) {
      if (verbose)
	printf ("GPTLsetutr: underlying wallclock timer = %s\n", 
		funclist[i].name);
      funcidx = i;
      return 0;
    }
  }
  return GPTLerror ("GPTLsetutr: unknown option %d\n", option);
}

/*
** GPTLinitialize (): Initialization routine must be called from single-threaded
**   region before any other timing routines may be called.  The need for this
**   routine could be eliminated if not targetting timing library for threaded
**   capability. 
**
** return value: 0 (success) or GPTLerror (failure)
*/

int GPTLinitialize (void)
{
  int i;          /* loop index */
  int t;          /* thread index */

  double t1, t2;  /* returned from underlying timer */

  if (initialized)
    return GPTLerror ("GPTLinitialize: has already been called\n");

  if (threadinit (&GPTLnthreads, &maxthreads) < 0)
    return GPTLerror ("GPTLinitialize: bad return from threadinit\n");

  if ((ticks_per_sec = sysconf (_SC_CLK_TCK)) == -1)
    return GPTLerror ("GPTLinitialize: sysconf (_SC_CLK_TCK) failed\n");

  /* Allocate space for global arrays */

  callstack     = (Timer ***)    GPTLallocate (maxthreads * sizeof (Timer **));
  stackidx      = (Nofalse *)    GPTLallocate (maxthreads * sizeof (Nofalse));
  timers        = (Timer **)     GPTLallocate (maxthreads * sizeof (Timer *));
  last          = (Timer **)     GPTLallocate (maxthreads * sizeof (Timer *));
  max_depth     = (int *)        GPTLallocate (maxthreads * sizeof (int));
  max_name_len  = (int *)        GPTLallocate (maxthreads * sizeof (int));
  hashtable     = (Hashentry **) GPTLallocate (maxthreads * sizeof (Hashentry *));

  /* Initialize array values */

  for (t = 0; t < maxthreads; t++) {
    max_depth[t]    = -1;
    max_name_len[t] = 0;
    callstack[t] = (Timer **) GPTLallocate (MAX_STACK * sizeof (Timer *));
    hashtable[t] = (Hashentry *) GPTLallocate (tablesize * sizeof (Hashentry));
    for (i = 0; i < tablesize; i++) {
      hashtable[t][i].nument = 0;
      hashtable[t][i].entries = 0;
    }

    /*
    ** Make a timer "GPTL_ROOT" to ensure no orphans, and to simplify printing.
    */

    timers[t] = (Timer *) GPTLallocate (sizeof (Timer));
    memset (timers[t], 0, sizeof (Timer));
    strcpy (timers[t]->name, "GPTL_ROOT");
    timers[t]->onflg = true;
    last[t] = timers[t];

    stackidx[t].val = 0;
    callstack[t][0] = timers[t];
    for (i = 1; i < MAX_STACK; i++)
      callstack[t][i] = 0;
  }

#ifdef HAVE_PAPI
  if (GPTL_PAPIinitialize (maxthreads, verbose, &nevents, eventlist) < 0)
    return GPTLerror ("GPTLinitialize: GPTL_PAPIinitialize failure\n");
#endif

  /* 
  ** Call init routine for underlying timing routine.
  */

  if ((*funclist[funcidx].funcinit)() < 0) {
    fprintf (stderr, "GPTLinitialize: failure initializing %s: "
	     "reverting underlying timer to %s\n", 
	     funclist[funcidx].name, funclist[0].name);
    funcidx = 0;
  }

  ptr2wtimefunc = funclist[funcidx].func;

  if (verbose) {
    t1 = (*ptr2wtimefunc) ();
    t2 = (*ptr2wtimefunc) ();
    if (t1 > t2)
      fprintf (stderr, "GPTLinitialize: negative delta-t=%g\n", t2-t1);
    printf ("Per call overhead est. t2-t1=%g should be near zero\n", t2-t1);
    printf ("Underlying wallclock timing routine is %s\n", funclist[funcidx].name);
  }

  initialized = true;
  return 0;
}

/*
** GPTLfinalize (): Finalization routine must be called from single-threaded
**   region. Free all malloc'd space
**
** return value: 0 (success) or GPTLerror (failure)
*/

int GPTLfinalize (void)
{
  int t;                /* thread index */
  Timer *ptr, *ptrnext; /* ll indices */

  if ( ! initialized)
    return GPTLerror ("GPTLfinalize: initialization was not completed\n");

  if (get_thread_num (&GPTLnthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLfinalize: must only be called by master thread\n");

  for (t = 0; t < maxthreads; ++t) {
    free (hashtable[t]);
    hashtable[t] = NULL;
    free (callstack[t]);
    for (ptr = timers[t]; ptr; ptr = ptrnext) {
      ptrnext = ptr->next;
      free (ptr);
    }
  }

  free (timers);
  free (max_depth);
  free (max_name_len);
  free (hashtable);
  free (callstack);
  free (stackidx);

  threadfinalize ();

#ifdef HAVE_PAPI
  GPTL_PAPIfinalize (maxthreads);
#endif

  /* Reset initial values set in GPTLinitialize */

  timers = 0;
  last = 0;
  GPTLnthreads = -1;
  maxthreads = -1;
  initialized = false;
  ref_gettimeofday = -1;
  ref_clock_gettime = -1;
  ref_papitime = -1;

  return 0;
}

/*
** GPTLstart_instr: start a timer (auto-instrumented)
**
** Input arguments:
**   self: function address
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLstart_instr (void *self)
{
  Timer *ptr;      /* linked list pointer */
  int t;           /* thread index (of this thread) */
  int indx;        /* hash table index */

  if (disabled)
    return 0;

  if ( ! initialized)
    return GPTLerror ("GPTLstart_instr self=%p: GPTLinitialize has not been called\n", self);

  if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
    return GPTLerror ("GPTLstart_instr: bad return from get_thread_num\n");

  /*
  ** If current depth exceeds a user-specified limit for print, just
  ** increment and return
  */

  if (stackidx[t].val >= depthlimit) {
    ++stackidx[t].val;
    return 0;
  }

  ptr = getentry_instr (hashtable[t], self, &indx);

  /* 
  ** Recursion => increment depth in recursion and return.  We need to return 
  ** because we don't want to restart the timer.  We want the reported time for
  ** the timer to reflect the outermost layer of recursion.
  */
  
  if (ptr && ptr->onflg) {
    ++ptr->recurselvl;
    return 0;
  }

  /*
  ** Increment stackidx[t] unconditionally. This is necessary to ensure the correct
  ** behavior when GPTLstop_instr decrements stackidx[t] unconditionally.
  */

  if (++stackidx[t].val > MAX_STACK-1)
    return GPTLerror ("GPTLstart_instr: stack too big\n");

  if ( ! ptr) {     /* Add a new entry and initialize */
    ptr = (Timer *) GPTLallocate (sizeof (Timer));
    memset (ptr, 0, sizeof (Timer));

    /*
    ** Need to save the address string for later conversion back to a real
    ** name by an offline tool.
    */

    snprintf (ptr->name, MAX_CHARS+1, "%lx", (unsigned long) self);
    ptr->address = self;

    if (update_ll_hash (ptr, t, indx) != 0)
      return GPTLerror ("GPTLstart_instr: update_ll_hash error\n");
  }

  if (update_parent_info (ptr, callstack[t], stackidx[t].val) != 0)
    return GPTLerror ("GPTLstart_instr: find_parent error\n");

  if (update_ptr (ptr, t) != 0)
    return GPTLerror ("GPTLstart_instr: update_ptr error\n");

  return (0);
}  

/*
** GPTLstart: start a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLstart (const char *name)               /* timer name */
{
  Timer *ptr;      /* linked list pointer */
  int t;           /* thread index (of this thread) */
  int indx;        /* hash table index */

#ifdef UNICOSMP
#ifndef SSP
  if (__streaming() == 0) return 0;  /* timers don't work in this situation so disable */
#endif
#endif

  if (disabled)
    return 0;

  if ( ! initialized)
    return GPTLerror ("GPTLstart name=%s: GPTLinitialize has not been called\n", name);

  if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
    return GPTLerror ("GPTLstart: bad return from get_thread_num\n");

  /*
  ** If current depth exceeds a user-specified limit for print, just
  ** increment and return
  */

  if (stackidx[t].val >= depthlimit) {
    ++stackidx[t].val;
    return 0;
  }

  /* 
  ** ptr will point to the requested timer in the current list,
  ** or NULL if this is a new entry 
  */

  ptr = getentry (hashtable[t], name, &indx);

  /* 
  ** Recursion => increment depth in recursion and return.  We need to return 
  ** because we don't want to restart the timer.  We want the reported time for
  ** the timer to reflect the outermost layer of recursion.
  */

  if (ptr && ptr->onflg) {
    ++ptr->recurselvl;
    return 0;
  }

  /*
  ** Increment stackidx[t] unconditionally. This is necessary to ensure the correct
  ** behavior when GPTLstop decrements stackidx[t] unconditionally.
  */

  if (++stackidx[t].val > MAX_STACK-1)
    return GPTLerror ("GPTLstart: stack too big\n");

  if ( ! ptr) { /* Add a new entry and initialize */
    ptr = (Timer *) GPTLallocate (sizeof (Timer));
    memset (ptr, 0, sizeof (Timer));

    strncpy (ptr->name, name, MAX_CHARS);
    ptr->name[MAX_CHARS] = '\0';

    if (update_ll_hash (ptr, t, indx) != 0)
      return GPTLerror ("GPTLstart: update_ll_hash error\n");
  }

  if (update_parent_info (ptr, callstack[t], stackidx[t].val) != 0)
    return GPTLerror ("GPTLstart: find_parent error\n");

  if (update_ptr (ptr, t) != 0)
    return GPTLerror ("GPTLstart: update_ptr error\n");

  return (0);
}

/*
** update_ll_hash: Update linked list and hash table.
**                 Called by GPTLstart and GPTLstart_instr
**
** Input arguments:
**   ptr:  pointer to timer
**   t:    thread index
**   indx: hash index
**
** Return value: 0 (success) or GPTLerror (failure)
*/

static inline int update_ll_hash (Timer *ptr, const int t, const int indx)
{
  int nchars;      /* number of chars */
  int nument;      /* number of entries */
  Timer **eptr;    /* for realloc */

  nchars = strlen (ptr->name);
  if (nchars > max_name_len[t])
    max_name_len[t] = nchars;

  last[t]->next = ptr;
  last[t] = ptr;
  ++hashtable[t][indx].nument;
  nument = hashtable[t][indx].nument;
  
  eptr = (Timer **) realloc (hashtable[t][indx].entries, nument * sizeof (Timer *));
  if ( ! eptr)
    return GPTLerror ("update_ll_hash: realloc error\n");

  hashtable[t][indx].entries           = eptr;
  hashtable[t][indx].entries[nument-1] = ptr;

  return 0;
}

/*
** update_ptr: Update timer contents. Called by GPTLstart and GPTLstart_instr
**
** Input arguments:
**   ptr:  pointer to timer
**   t:    thread index
**
** Return value: 0 (success) or GPTLerror (failure)
*/

static inline int update_ptr (Timer *ptr, const int t)
{
  double tp2;    /* time stamp */

  ptr->onflg = true;

  if (cpustats.enabled && get_cpustamp (&ptr->cpu.last_utime, &ptr->cpu.last_stime) < 0)
    return GPTLerror ("update_ptr: get_cpustamp error");
  
  if (wallstats.enabled) {
    tp2 = (*ptr2wtimefunc) ();
    ptr->wall.last = tp2;
  }

#ifdef HAVE_PAPI
  if (dousepapi && GPTL_PAPIstart (t, &ptr->aux) < 0)
    return GPTLerror ("update_ptr: error from GPTL_PAPIstart\n");
#endif
  return 0;
}

/*
** update_parent_info: update info about parent, and in the parent about this child
**
** Input arguments:
**   ptr:  pointer to timer
**   callstackt: callstack for this thread
**   stackidxt:  stack index for this thread
**
** Return value: 0 (success) or GPTLerror (failure)
*/

static inline int update_parent_info (Timer *ptr, Timer **callstackt, int stackidxt) 
{
  int n;             /* loop index through known parents */
  Timer *pptr;       /* pointer to parent */
  Timer **pptrtmp;   /* for realloc parent pointer array */
  int nparent;       /* number of parents */
  int *parent_count; /* number of times parent invoked this child */

  if ( ! ptr )
    return -1;

  if (stackidxt < 0)
    return GPTLerror ("update_parent_info: called with negative stackidx\n");

  callstackt[stackidxt] = ptr;

  /* 
  ** If the region has no parent, bump its orphan count (should never happen since "GPTL_ROOT" added).
  ** Don't set depth=0, since the region may have had a parent earlier.
  ** If the region is always an orphan, the initial memset already set its
  ** depth to 0.
  */

  if (stackidxt == 0) {
    ++ptr->norphan;
    return 0;
  }

  pptr = callstackt[stackidxt-1];

  /* If this parent occurred before, bump its count */

  for (n = 0; n < ptr->nparent; ++n) {
    if (ptr->parent[n] == pptr) {
      ++ptr->parent_count[n];
      break;
    }
  }

  /* If this is a new parent, update info */

  if (n == ptr->nparent) {
    ++ptr->nparent;
    nparent = ptr->nparent;
    pptrtmp = (Timer **) realloc (ptr->parent, nparent * sizeof (Timer *));
    if ( ! pptrtmp)
      return GPTLerror ("update_parent_info: realloc error pptrtmp nparent=%d\n", nparent);

    ptr->parent = pptrtmp;
    ptr->parent[nparent-1] = pptr;
    parent_count = (int *) realloc (ptr->parent_count, nparent * sizeof (int));
    if ( ! parent_count)
      return GPTLerror ("update_parent_info: realloc error parent_count nparent=%d\n", nparent);

    ptr->parent_count = parent_count;
    ptr->parent_count[nparent-1] = 1;
  }

  return 0;
}

/*
** GPTLstop_instr: stop a timer (auto-instrumented)
**
** Input arguments:
**   self: function address
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLstop_instr (void *self)
{
  double tp1 = 0.0;          /* time stamp */
  Timer *ptr;                /* linked list pointer */

  int t;                     /* thread number for this process */
  int indx;                  /* index into hash table */

  long usr = 0;              /* user time (returned from get_cpustamp) */
  long sys = 0;              /* system time (returned from get_cpustamp) */

#ifdef UNICOSMP
#ifndef SSP
  if (__streaming() == 0) return 0;  /* timers don't work in this situation so disable */
#endif
#endif

  if (disabled)
    return 0;

  if ( ! initialized)
    return GPTLerror ("GPTLstop_instr: GPTLinitialize has not been called\n");

  /* Get the timestamp */
    
  if (wallstats.enabled) {
    tp1 = (*ptr2wtimefunc) ();
  }

  if (cpustats.enabled && get_cpustamp (&usr, &sys) < 0)
    return GPTLerror (0);

  if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
    return GPTLerror ("GPTLstop_instr\n");

  /*
  ** If current depth exceeds a user-specified limit for print, just
  ** decrement and return
  */

  if (stackidx[t].val > depthlimit) {
    --stackidx[t].val;
    return 0;
  }

  ptr = getentry_instr (hashtable[t], self, &indx);

  if ( ! ptr) 
    return GPTLerror ("GPTLstop_instr: timer for %p had not been started.\n", self);

  if ( ! ptr->onflg )
    return GPTLerror ("GPTLstop_instr: timer %s was already off.\n",ptr->name);

  ++ptr->count;

  /* 
  ** Recursion => decrement depth in recursion and return.  We need to return
  ** because we don't want to stop the timer.  We want the reported time for
  ** the timer to reflect the outermost layer of recursion.
  */

  if (ptr->recurselvl > 0) {
    ++ptr->nrecurse;
    --ptr->recurselvl;
    return 0;
  }

  if (update_stats (ptr, tp1, usr, sys, t) != 0)
    return GPTLerror ("GPTLstop_instr: update_stats\n");

  return 0;
}

/*
** GPTLstop: stop a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or -1 (failure)
*/

int GPTLstop (const char *name)               /* timer name */
{
  double tp1 = 0.0;          /* time stamp */
  Timer *ptr;                /* linked list pointer */
  int t;                     /* thread number for this process */
  int indx;                  /* index into hash table */
  long usr = 0;              /* user time (returned from get_cpustamp) */
  long sys = 0;              /* system time (returned from get_cpustamp) */

#ifdef UNICOSMP
#ifndef SSP
  if (__streaming() == 0) return 0;  /* timers don't work in this situation so disable */
#endif
#endif

  if (disabled)
    return 0;

  if ( ! initialized)
    return GPTLerror ("GPTLstop: GPTLinitialize has not been called\n");

  /* Get the timestamp */
    
  if (wallstats.enabled) {
    tp1 = (*ptr2wtimefunc) ();
  }

  if (cpustats.enabled && get_cpustamp (&usr, &sys) < 0)
    return GPTLerror (0);

  if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
    return GPTLerror ("GPTLstop\n");

  /*
  ** If current depth exceeds a user-specified limit for print, just
  ** decrement and return
  */

  if (stackidx[t].val > depthlimit) {
    --stackidx[t].val;
    return 0;
  }

  if ( ! (ptr = getentry (hashtable[t], name, &indx)))
    return GPTLerror ("GPTLstop: timer for %s had not been started.\n", name);

  if ( ! ptr->onflg )
    return GPTLerror ("GPTLstop: timer %s was already off.\n",ptr->name);

  ++ptr->count;

  /* 
  ** Recursion => decrement depth in recursion and return.  We need to return
  ** because we don't want to stop the timer.  We want the reported time for
  ** the timer to reflect the outermost layer of recursion.
  */

  if (ptr->recurselvl > 0) {
    ++ptr->nrecurse;
    --ptr->recurselvl;
    return 0;
  }

  if (update_stats (ptr, tp1, usr, sys, t) != 0)
    return GPTLerror ("GPTLstop: update_stats\n");

  return 0;
}

/*
** update_stats: update stats inside ptr. Called by GPTLstop, GPTLstop_instr
**
** Input arguments:
**   ptr: pointer to timer
**   tp1: input time stapm
**   usr: user time
**   sys: system time
**   t: thread index
**
** Return value: 0 (success) or GPTLerror (failure)
*/

static inline int update_stats (Timer *ptr, 
				const double tp1, 
				const long usr, 
				const long sys,
				const int t)
{
  double delta;       /* difference */

  ptr->onflg = false;
  --stackidx[t].val;
  if (stackidx[t].val < -1) {
    stackidx[t].val = -1;
    return GPTLerror ("update_stats: tree depth has become negative.\n");
  }

#ifdef HAVE_PAPI
  if (dousepapi && GPTL_PAPIstop (t, &ptr->aux) < 0)
    return GPTLerror ("update_stats: error from GPTL_PAPIstop\n");
#endif

  if (wallstats.enabled) {
    delta = tp1 - ptr->wall.last;
    ptr->wall.accum += delta;

    if (delta < 0.) {
      fprintf (stderr, "update_stats: negative delta=%g\n", delta);
    }

    if (ptr->count == 1) {
      ptr->wall.max = delta;
      ptr->wall.min = delta;
    } else {
      if (delta > ptr->wall.max)
	ptr->wall.max = delta;
      if (delta < ptr->wall.min)
	ptr->wall.min = delta;
    }
  }

  if (cpustats.enabled) {
    ptr->cpu.accum_utime += usr - ptr->cpu.last_utime;
    ptr->cpu.accum_stime += sys - ptr->cpu.last_stime;
    ptr->cpu.last_utime   = usr;
    ptr->cpu.last_stime   = sys;
  }
  return 0;
}

/*
** GPTLenable: enable timers
**
** Return value: 0 (success)
*/

int GPTLenable (void)
{
  disabled = false;
  return (0);
}

/*
** GPTLdisable: disable timers
**
** Return value: 0 (success)
*/

int GPTLdisable (void)
{
  disabled = true;
  return (0);
}

/*
** GPTLstamp: Compute timestamp of usr, sys, and wallclock time (seconds)
**
** Output arguments:
**   wall: wallclock
**   usr:  user time
**   sys:  system time
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLstamp (double *wall, double *usr, double *sys)
{
  struct tms buf;            /* argument to times */

  if ( ! initialized)
    return GPTLerror ("GPTLstamp: GPTLinitialize has not been called\n");

#ifdef HAVE_TIMES
  *usr = 0;
  *sys = 0;

  if (times (&buf) == -1)
    return GPTLerror ("GPTLstamp: times() failed. Results bogus\n");

  *usr = buf.tms_utime / (double) ticks_per_sec;
  *sys = buf.tms_stime / (double) ticks_per_sec;
#endif
  *wall = (*ptr2wtimefunc) ();
  return 0;
}

/*
** GPTLreset: reset all known timers to 0
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLreset (void)
{
  int t;             /* index over threads */
  Timer *ptr;        /* linked list index */

  if ( ! initialized)
    return GPTLerror ("GPTLreset: GPTLinitialize has not been called\n");

  /* Only allow the master thread to reset timers */

  if (get_thread_num (&GPTLnthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLreset: must only be called by master thread\n");

  for (t = 0; t < GPTLnthreads; t++) {
    for (ptr = timers[t]; ptr; ptr = ptr->next) {
      ptr->onflg = false;
      ptr->count = 0;
      memset (&ptr->wall, 0, sizeof (ptr->wall));
      memset (&ptr->cpu, 0, sizeof (ptr->cpu));
#ifdef HAVE_PAPI
      memset (&ptr->aux, 0, sizeof (ptr->aux));
#endif
    }
  }

  if (verbose)
    printf ("GPTLreset: accumulators for all timers set to zero\n");

  return 0;
}

/* 
** GPTLpr: Print values of all timers
**
** Input arguments:
**   id: integer to append to string "timing."
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLpr (const int id)   /* output file will be named "timing.<id>" */
{
  char outfile[14];         /* name of output file: timing.xxxxxx */

  if (id < 0 || id > 999999)
    return GPTLerror ("GPTLpr: bad id=%d for output file. Must be >= 0 and < 1000000\n",
		      id);

  sprintf (outfile, "timing.%d", id);

  if (GPTLpr_file (outfile) != 0)
    return GPTLerror ("GPTLpr: Error in GPTLpr_file\n");

  return 0;
}

/* 
** GPTLpr_file: Print values of all timers
**
** Input arguments:
**   outfile: Name of output file to write
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLpr_file (const char *outfile) /* output file to write */
{
  FILE *fp;                 /* file handle to write to */
  Timer *ptr;               /* walk through master thread linked list */
  Timer *tptr;              /* walk through slave threads linked lists */
  Timer sumstats;           /* sum of same timer stats over threads */
  int i, ii, n, t;          /* indices */
  int totent;               /* per-thread collision count (diagnostic) */
  int nument;               /* per-index collision count (diagnostic) */
  int totlen;               /* length for malloc */
  unsigned long totcount;   /* total timer invocations */
  unsigned long totrecurse; /* total recursive timer invocations */
  char *outpath;            /* path to output file: outdir/timing.xxxxxx */
  float *sum;               /* sum of overhead values (per thread) */
  float osum;               /* sum of overhead over threads */
  double utr_overhead;      /* overhead of calling underlying timing routine */
  double tot_overhead;      /* utr_overhead + papi overhead */
  double papi_overhead = 0; /* overhead of reading papi counters */
  bool found;               /* jump out of loop when name found */
  bool foundany;            /* whether summation print necessary */
  bool first;               /* flag 1st time entry found */
  /*
  ** Diagnostics for collisions and GPTL memory usage
  */
  int num_zero;             /* number of buckets with 0 collisions */
  int num_one;              /* number of buckets with 1 collision */
  int num_two;              /* number of buckets with 2 collisions */
  int num_more;             /* number of buckets with more than 2 collisions */
  int most;                 /* biggest collision count */
  int numtimers;            /* number of timers */
  int hashmem;              /* hash table memory usage */
  int regionmem;            /* timer memory usage */
  int papimem;              /* PAPI stats memory usage */
  int pchmem;               /* parent/child array memory usage */
  int gptlmem;              /* total GPTL memory usage estimate */

  if ( ! initialized)
    return GPTLerror ("GPTLpr_file: GPTLinitialize() has not been called\n");

  if (get_thread_num (&GPTLnthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLpr_file: must only be called by master thread\n");

  /* 2 is for "/" plus null */
  if (outdir)
    totlen = strlen (outdir) + strlen (outfile) + 2; 
  else
    totlen = strlen (outfile) + 2; 

  outpath = (char *) GPTLallocate (totlen);

  if (outdir) {
     strcpy (outpath, outdir);
     strcat (outpath, "/");
     strcat (outpath, outfile);
  } else {
     strcpy (outpath, outfile);
  }

  if ( ! (fp = fopen (outpath, "w")))
    fp = stderr;

  free (outpath);

  fprintf (fp, "$Id: gptl.c,v 1.153 2010/04/13 21:41:50 rosinski Exp $\n");

#ifdef HAVE_NANOTIME
  if (funcidx == GPTLnanotime)
    fprintf (fp, "Clock rate = %f MHz\n", cpumhz);
#endif

#ifdef HAVE_PAPI
  if (dousepapi) {
    if (GPTL_PAPIis_multiplexed ())
      fprintf (fp, "PAPI event multiplexing was ON\n");
    else
      fprintf (fp, "PAPI event multiplexing was OFF\n");
    GPTL_PAPIprintenabled (fp);
  }
#endif

  /*
  ** Estimate underlying timing routine overhead
  */

  utr_overhead = utr_getoverhead ();
  fprintf (fp, "Underlying timing routine was %s.\n", funclist[funcidx].name);
  fprintf (fp, "Per-call utr overhead est: %g sec.\n", utr_overhead);
#ifdef HAVE_PAPI
  if (dousepapi) {
    double t1, t2;
    t1 = (*ptr2wtimefunc) ();
    read_counters100 ();
    t2 = (*ptr2wtimefunc) ();
    papi_overhead = 0.01 * (t2 - t1);
    fprintf (fp, "Per-call PAPI overhead est: %g sec.\n", papi_overhead);
  }
#endif
  tot_overhead = utr_overhead + papi_overhead;
  if (dopr_preamble) {
    fprintf (fp, "If overhead stats are printed, roughly half the estimated number is\n");
    fprintf (fp, "embedded in the wallclock stats for each timer.\n");
    fprintf (fp, "Print method was %s.\n", methodstr (method));
#ifdef ENABLE_PMPI
    fprintf (fp, "If a AVG_MPI_BYTES field is present, it is an estimate of the per-call "
	     "average number of bytes handled by that process.\n");
    fprintf (fp, "If timers beginning with sync_ are present, it means MPI synchronization "
	     "was turned on.\n");
#endif
    fprintf (fp, "If a \'%%_of\' field is present, it is w.r.t. the first timer for thread 0.\n");
    fprintf (fp, "If a \'e6_per_sec\' field is present, it is in millions of PAPI counts per sec.\n\n");
    fprintf (fp, "A '*' in column 1 below means the timer had multiple parents, though the\n");
    fprintf (fp, "values printed are for all calls. Further down the listing is more detailed\n");
    fprintf (fp, "information about multiple parents. Look for 'Multiple parent info'\n\n");
  }

  sum = (float *) GPTLallocate (GPTLnthreads * sizeof (float));
  
  for (t = 0; t < GPTLnthreads; ++t) {

    /*
    ** Construct tree for printing timers in parent/child form. get_max_depth() must be called 
    ** AFTER construct_tree() because it relies on the per-parent children arrays being complete.
    */

    if (construct_tree (timers[t], method) != 0)
      printf ("GPTLpr_file: failure from construct_tree: output will be incomplete\n");
    max_depth[t] = get_max_depth (timers[t], 0);

    if (t > 0)
      fprintf (fp, "\n");
    fprintf (fp, "Stats for thread %d:\n", t);

    for (n = 0; n < max_depth[t]+1; ++n)    /* +1 to always indent timer name */
      fprintf (fp, "  ");
    for (n = 0; n < max_name_len[t]; ++n) /* longest timer name */
      fprintf (fp, " ");
    fprintf (fp, "Called  Recurse ");

    /* Print strings for enabled timer types */

    if (cpustats.enabled)
      fprintf (fp, "%s", cpustats.str);
    if (wallstats.enabled) {
      fprintf (fp, "%s", wallstats.str);
      if (percent && timers[0]->next)
	fprintf (fp, "%%_of_%5.5s ", timers[0]->next->name);
      if (overheadstats.enabled)
	fprintf (fp, "%s", overheadstats.str);
    }

#ifdef ENABLE_PMPI
    fprintf (fp, "AVG_MPI_BYTES ");
#endif

#ifdef HAVE_PAPI
    GPTL_PAPIprstr (fp);
#endif

    fprintf (fp, "\n");        /* Done with titles, now print stats */

    /*
    ** Print call tree and stats via recursive routine. "-1" is flag to
    ** avoid printing dummy outermost timer, and initialize the depth.
    */

    printself_andchildren (timers[t], fp, t, -1, tot_overhead);

    /* 
    ** Sum of overhead across timers is meaningful.
    ** Factor of 2 is because there are 2 utr calls per start/stop pair.
    */

    sum[t]     = 0;
    totcount   = 0;
    totrecurse = 0;
    for (ptr = timers[t]->next; ptr; ptr = ptr->next) {
      sum[t]     += ptr->count * 2 * tot_overhead;
      totcount   += ptr->count;
      totrecurse += ptr->nrecurse;
    }
    if (wallstats.enabled && overheadstats.enabled)
      fprintf (fp, "Overhead sum          = %9.3g wallclock seconds\n", sum[t]);
    if (totcount < PRTHRESH)
      fprintf (fp, "Total calls           = %lu\n", totcount);
    else
      fprintf (fp, "Total calls           = %8.1e\n", (float) totcount);
    fprintf (fp, "Total recursive calls = %lu\n", totrecurse);
  }

  /* Print per-name stats for all threads */

  if (dopr_threadsort && GPTLnthreads > 1) {
    fprintf (fp, "\nSame stats sorted by timer for threaded regions:\n");
    fprintf (fp, "Thd ");

    for (n = 0; n < max_name_len[0]; ++n) /* longest timer name */
      fprintf (fp, " ");

    fprintf (fp, "Called  Recurse ");

    if (cpustats.enabled)
      fprintf (fp, "%s", cpustats.str);
    if (wallstats.enabled) {
      fprintf (fp, "%s", wallstats.str);
      if (percent && timers[0]->next)
	fprintf (fp, "%%_of_%5.5s ", timers[0]->next->name);
      if (overheadstats.enabled)
	fprintf (fp, "%s", overheadstats.str);
    }

#ifdef HAVE_PAPI
    GPTL_PAPIprstr (fp);
#endif

    fprintf (fp, "\n");

    /* Start at next to skip dummy */

    for (ptr = timers[0]->next; ptr; ptr = ptr->next) {
      
      /* 
      ** To print sum stats, first create a new timer then copy thread 0
      ** stats into it. then sum using "add", and finally print.
      */

      foundany = false;
      first = true;
      sumstats = *ptr;
      for (t = 1; t < GPTLnthreads; ++t) {
	found = false;
	for (tptr = timers[t]->next; tptr && ! found; tptr = tptr->next) {
	  if (STRMATCH (ptr->name, tptr->name)) {

	    /* Only print thread 0 when this timer found for other threads */

	    if (first) {
	      first = false;
	      fprintf (fp, "%3.3d ", 0);
	      printstats (ptr, fp, 0, 0, false, tot_overhead);
	    }

	    found = true;
	    foundany = true;
	    fprintf (fp, "%3.3d ", t);
	    printstats (tptr, fp, 0, 0, false, tot_overhead);
	    add (&sumstats, tptr);
	  }
	}
      }

      if (foundany) {
	fprintf (fp, "SUM ");
	printstats (&sumstats, fp, 0, 0, false, tot_overhead);
	fprintf (fp, "\n");
      }
    }

    /* Repeat overhead print in loop over threads */

    if (wallstats.enabled && overheadstats.enabled) {
      osum = 0.;
      for (t = 0; t < GPTLnthreads; ++t) {
	fprintf (fp, "OVERHEAD.%3.3d (wallclock seconds) = %9.3g\n", t, sum[t]);
	osum += sum[t];
      }
      fprintf (fp, "OVERHEAD.SUM (wallclock seconds) = %9.3g\n", osum);
    }
  }

  /* Print info about timers with multiple parents */

  if (dopr_multparent) {
    for (t = 0; t < GPTLnthreads; ++t) {
      fprintf (fp, "\nMultiple parent info (if any) for thread %d:\n", t);
      if (dopr_preamble && t == 0) {
	fprintf (fp, "Columns are count and name for the listed child\n");
	fprintf (fp, "Rows are each parent, with their common child being the last entry, "
		 "which is indented\n");
	fprintf (fp, "Count next to each parent is the number of times it called the "
		 "child\n");
	fprintf (fp, "Count next to child is total number of times it was called by the "
		 "listed parents\n\n");
      }

      for (ptr = timers[t]->next; ptr; ptr = ptr->next)
	if (ptr->nparent > 1)
	  print_multparentinfo (fp, ptr);
    }
  }

  /* Print hash table stats */

  if (dopr_collision) {
    for (t = 0; t < GPTLnthreads; t++) {
      first = true;
      totent   = 0;
      num_zero = 0;
      num_one  = 0;
      num_two  = 0;
      num_more = 0;
      most     = 0;
      numtimers= 0;

      for (i = 0; i < tablesize; i++) {
	nument = hashtable[t][i].nument;
	if (nument > 1) {
	  totent += nument-1;
	  if (first) {
	    first = false;
	    fprintf (fp, "\nthread %d had some hash collisions:\n", t);
	  }
	  fprintf (fp, "hashtable[%d][%d] had %d entries:", t, i, nument);
	  for (ii = 0; ii < nument; ii++)
	    fprintf (fp, " %s", hashtable[t][i].entries[ii]->name);
	  fprintf (fp, "\n");
	}
	switch (nument) {
	case 0:
	  ++num_zero;
	  break;
	case 1:
	  ++num_one;
	  break;
	case 2:
	  ++num_two;
	  break;
	default:
	  ++num_more;
	  break;
	}
	most = MAX (most, nument);
	numtimers += nument;
      }

      if (totent > 0) {
	fprintf (fp, "Total collisions thread %d = %d\n", t, totent);
	fprintf (fp, "Entry information:\n");
	fprintf (fp, "num_zero = %d num_one = %d num_two = %d num_more = %d\n",
		 num_zero, num_one, num_two, num_more);
	fprintf (fp, "Most = %d\n", most);
      }

      /* Stats on GPTL memory usage */

      hashmem = sizeof (Hashentry) * tablesize;
      regionmem = numtimers * sizeof (Timer);
#ifdef HAVE_PAPI
      papimem = numtimers * sizeof (Papistats);
#else
      papimem = 0.;
#endif
      pchmem = 0.;
      for (ptr = timers[t]->next; ptr; ptr = ptr->next)
	pchmem += (sizeof (Timer *)) * (ptr->nchildren + ptr->nparent);

      gptlmem = hashmem + regionmem + pchmem;
      fprintf (fp, "Thread %d total memory usage=%g KB\n", t, gptlmem*1.e-3);
      fprintf (fp, "hashmem = %g KB\n" 
	       "regionmem = %g KB (papimem portion = %g KB)\n"
	       "parent/child arrays = %g KB\n",
	       hashmem*.001, regionmem*.001, papimem*.001, pchmem*.001);
    }
  }

  print_threadmapping (GPTLnthreads, fp);

  free (sum);

  if (fclose (fp) != 0)
    fprintf (stderr, "Attempt to close %s failed\n", outfile);

  pr_has_been_called = true;
  return 0;
}

/* 
** construct_tree: Build the parent->children tree starting with knowledge of
**                 parent list for each child.
**
** Input arguments:
**   timerst: Linked list of timers
**   method:  method to be used to define the links
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int construct_tree (Timer *timerst, Method method)
{
  Timer *ptr;       /* loop through linked list */
  Timer *pptr = 0;  /* parent (init to NULL to avoid compiler warning) */
  int nparent;      /* number of parents */
  int maxcount;     /* max calls by a single parent */
  int n;            /* loop over nparent */

  /*
  ** Walk the linked list to build the parent-child tree, using whichever
  ** mechanism is in place. newchild() will prevent loops.
  */

  for (ptr = timerst; ptr; ptr = ptr->next) {
    switch (method) {
    case GPTLfirst_parent:
      if (ptr->nparent > 0) {
	pptr = ptr->parent[0];
	if (newchild (pptr, ptr) != 0);
      }
      break;
    case GPTLlast_parent:
      if (ptr->nparent > 0) {
	nparent = ptr->nparent;
	pptr = ptr->parent[nparent-1];
	if (newchild (pptr, ptr) != 0);
      }
      break;
    case GPTLmost_frequent:
      maxcount = 0;
      for (n = 0; n < ptr->nparent; ++n) {
	if (ptr->parent_count[n] > maxcount) {
	  pptr = ptr->parent[n];
	  maxcount = ptr->parent_count[n];
	}
      }
      if (maxcount > 0) {   /* not an orphan */
	if (newchild (pptr, ptr) != 0);
      }
      break;
    case GPTLfull_tree:
      /* 
      ** Careful: this one can create *lots* of output!
      */
      for (n = 0; n < ptr->nparent; ++n) {
	pptr = ptr->parent[n];
	if (newchild (pptr, ptr) != 0);
      }
      break;
    default:
      return GPTLerror ("construct_tree: method %d is not known\n", method);
    }
  }
  return 0;
}

/* 
** methodstr: Return a pointer to a string which represents the method
**
** Input arguments:
**   method: method type
*/

static char *methodstr (Method method)
{
  if (method == GPTLfirst_parent)
    return "first_parent";
  else if (method == GPTLlast_parent)
    return "last_parent";
  else if (method == GPTLmost_frequent)
    return "most_frequent";
  else if (method == GPTLfull_tree)
    return "full_tree";
  else
    return "Unknown";
}

/* 
** newchild: Add an entry to the children list of parent. Use function
**   is_descendant() to prevent infinite loops. 
**
** Input arguments:
**   parent: parent node
**   child:  child to be added
**
** Return value: 0 (success) or GPTLerror (failure)
*/

static int newchild (Timer *parent, Timer *child)
{
  int nchildren;     /* number of children (temporary) */
  Timer **chptr;     /* array of pointers to children */

  if (parent == child)
    return GPTLerror ("newchild: child %s can't be a parent of itself\n", child->name);

  /*
  ** To guarantee no loops, ensure that proposed parent isn't already a descendant of 
  ** proposed child
  */

  if (is_descendant (child, parent)) {
    return GPTLerror ("newchild: loop detected: NOT adding %s to descendant list of %s. "
		      "Proposed parent is in child's descendant path.\n",
		      child->name, parent->name);
  }

  /* Safe to add the child to the parent's list of children */

  ++parent->nchildren;
  nchildren = parent->nchildren;
  chptr = (Timer **) realloc (parent->children, nchildren * sizeof (Timer *));
  if ( ! chptr)
    return GPTLerror ("newchild: realloc error\n");
  parent->children = chptr;
  parent->children[nchildren-1] = child;

  return 0;
}

/* 
** get_max_depth: Determine the maximum call tree depth by traversing the
**   tree recursively
**
** Input arguments:
**   ptr:        Starting timer
**   startdepth: current depth when function invoked 
**
** Return value: maximum depth
*/

static int get_max_depth (const Timer *ptr, const int startdepth)
{
  int maxdepth = startdepth;
  int depth;
  int n;

  for (n = 0; n < ptr->nchildren; ++n)
    if ((depth = get_max_depth (ptr->children[n], startdepth+1)) > maxdepth)
      maxdepth = depth;

  return maxdepth;
}

/* 
** num_descendants: Determine the number of descendants of a timer by traversing
**   the tree recursively. This function is not currently used. It could be
**   useful in a pruning algorithm
**
** Input arguments:
**   ptr: Starting timer
**
** Return value: number of descendants
*/

static int num_descendants (Timer *ptr)
{
  int n;

  ptr->num_desc = ptr->nchildren;
  for (n = 0; n < ptr->nchildren; ++n) {
    ptr->num_desc += num_descendants (ptr->children[n]);
  }
  return ptr->num_desc;
}

/* 
** is_descendant: Determine whether node2 is in the descendant list for
**   node1
**
** Input arguments:
**   node1: starting node for recursive search
**   node2: node to be searched for
**
** Return value: true or false
*/

static int is_descendant (const Timer *node1, const Timer *node2)
{
  int n;

  /* Breadth before depth for efficiency */

  for (n = 0; n < node1->nchildren; ++n)
    if (node1->children[n] == node2)
      return 1;

  for (n = 0; n < node1->nchildren; ++n)
    if (is_descendant (node1->children[n], node2))
      return 1;

  return 0;
}

/* 
** printstats: print a single timer
**
** Input arguments:
**   timer:        timer for which to print stats
**   fp:           file descriptor to write to
**   t:            thread number
**   depth:        depth to indent timer
**   doindent:     whether indenting will be done
**   tot_overhead: underlying timing routine overhead
*/

static void printstats (const Timer *timer,
			FILE *fp,
			const int t,
			const int depth,
			const bool doindent,
			const double tot_overhead)
{
  int i;               /* index */
  int indent;          /* index for indenting */
  int extraspace;      /* for padding to length of longest name */
  float fusr;          /* user time as float */
  float fsys;          /* system time as float */
  float usrsys;        /* usr + sys */
  float elapse;        /* elapsed time */
  float wallmax;       /* max wall time */
  float wallmin;       /* min wall time */
  float ratio;         /* percentage calc */

  if (timer->onflg && verbose)
    fprintf (stderr, "printstats: timer %s had not been turned off\n", timer->name);

  /* Flag regions having multiple parents with a "*" in column 1 */

  if (doindent) {
    if (timer->nparent > 1)
      fprintf (fp, "* ");
    else
      fprintf (fp, "  ");

    /* Indent to depth of this timer */

    for (indent = 0; indent < depth; ++indent)
      fprintf (fp, "  ");
  }

  fprintf (fp, "%s", timer->name);

  /* Pad to length of longest name */

  extraspace = max_name_len[t] - strlen (timer->name);
  for (i = 0; i < extraspace; ++i)
    fprintf (fp, " ");

  /* Pad to max indent level */

  if (doindent)
    for (indent = depth; indent < max_depth[t]; ++indent)
      fprintf (fp, "  ");

  if (timer->count < PRTHRESH) {
    if (timer->nrecurse > 0)
      fprintf (fp, "%8lu %6lu ", timer->count, timer->nrecurse);
    else
      fprintf (fp, "%8lu    -   ", timer->count);
  } else {
    if (timer->nrecurse > 0)
      fprintf (fp, "%8.1e %6.0e ", (float) timer->count, (float) timer->nrecurse);
    else
      fprintf (fp, "%8.1e    -   ", (float) timer->count);
  }

  if (cpustats.enabled) {
    fusr = timer->cpu.accum_utime / (float) ticks_per_sec;
    fsys = timer->cpu.accum_stime / (float) ticks_per_sec;
    usrsys = fusr + fsys;
    fprintf (fp, "%9.3f %9.3f %9.3f ", fusr, fsys, usrsys);
  }

  if (wallstats.enabled) {
    elapse = timer->wall.accum;
    wallmax = timer->wall.max;
    wallmin = timer->wall.min;
    fprintf (fp, "%9.3f %9.3f %9.3f ", elapse, wallmax, wallmin);

    if (percent && timers[0]->next) {
      ratio = 0.;
      if (timers[0]->next->wall.accum > 0.)
	ratio = (timer->wall.accum * 100.) / timers[0]->next->wall.accum;
      fprintf (fp, " %9.2f ", ratio);
    }

    /*
    ** Factor of 2 is because there are 2 utr calls per start/stop pair.
    */

    if (overheadstats.enabled) {
      fprintf (fp, "%13.3f ", timer->count * 2 * tot_overhead);
    }
  }

#ifdef ENABLE_PMPI
  if (timer->nbytes == 0.)
    fprintf (fp, "      -       ");
  else
    fprintf (fp, "%13.3e ", timer->nbytes / timer->count);
#endif
  
#ifdef HAVE_PAPI
  GPTL_PAPIpr (fp, &timer->aux, t, timer->count, timer->wall.accum);
#endif

  fprintf (fp, "\n");
}

/* 
** print_multparentinfo: 
**
** Input arguments:
** Input/output arguments:
*/
void print_multparentinfo (FILE *fp, 
			   Timer *ptr)
{
  int n;

  if (ptr->norphan > 0) {
    if (ptr->norphan < PRTHRESH)
      fprintf (fp, "%8u %-32s\n", ptr->norphan, "ORPHAN");
    else
      fprintf (fp, "%8.1e %-32s\n", (float) ptr->norphan, "ORPHAN");
  }

  for (n = 0; n < ptr->nparent; ++n) {
    if (ptr->parent_count[n] < PRTHRESH)
      fprintf (fp, "%8d %-32s\n", ptr->parent_count[n], ptr->parent[n]->name);
    else
      fprintf (fp, "%8.1e %-32s\n", (float) ptr->parent_count[n], ptr->parent[n]->name);
  }

  if (ptr->count < PRTHRESH)
    fprintf (fp, "%8lu   %-32s\n\n", ptr->count, ptr->name);
  else
    fprintf (fp, "%8.1e   %-32s\n\n", (float) ptr->count, ptr->name);
}

/* 
** add: add the contents of tin to tout
**
** Input arguments:
**   tin:  input timer
** Input/output arguments:
**   tout: output timer summed into
*/

static void add (Timer *tout,   
		 const Timer *tin)
{
  if (wallstats.enabled) {
    tout->count      += tin->count;
    tout->wall.accum += tin->wall.accum;
    
    tout->wall.max = MAX (tout->wall.max, tin->wall.max);
    tout->wall.min = MIN (tout->wall.min, tin->wall.min);
  }

  if (cpustats.enabled) {
    tout->cpu.accum_utime += tin->cpu.accum_utime;
    tout->cpu.accum_stime += tin->cpu.accum_stime;
  }
#ifdef HAVE_PAPI
  GPTL_PAPIadd (&tout->aux, &tin->aux);
#endif
}

#ifdef HAVE_MPI
/* 
** GPTLpr_summary: When MPI enabled, gather and print summary stats across 
**                 threads and MPI tasks
**
** Input arguments:
**   comm: commuicator (e.g. MPI_COMM_WORLD). If zero, use MPI_COMM_WORLD
*/

int GPTLpr_summary (MPI_Comm comm)
{
  const char *outfile = "timing.summary";
  int iam = 0;                     /* MPI rank: default master */
  int n;                           /* index */
  int extraspace;                  /* for padding to length of longest name */
  Summarystats summarystats;       /* stats to be printed */
  Timer *ptr;                      /* timer */
  FILE *fp = 0;                    /* output file */

  int nproc;                                /* number of procs in MPI communicator */
  int p;                                    /* process index */
  int ret;                                  /* return code */
  const int nbytes = sizeof (Summarystats); /* number of bytes to be sent */
  const int tag = 99;                       /* tag for MPI message */
  char name[MAX_CHARS+1];                   /* timer name requested by master */
  Summarystats summarystats_slave;          /* stats sent to master */
  MPI_Status status;                        /* required by MPI_Recv */
  char emptystring[MAX_CHARS+1];            /* needs to be as big as name */

  /*
  ** Ensure that emptystring starts with \0
  */

  emptystring[0] = '\0';

  if (((int) comm) == 0)
    comm = MPI_COMM_WORLD;

  ret = MPI_Comm_rank (comm, &iam);
  ret = MPI_Comm_size (comm, &nproc);

  /*
  ** Master loops over thread 0 timers. Each process gathers stats for its threads,
  ** master gathers these results.
  ** End of requests signaled by requesting the null string.
  */

  if (iam == 0) {
    if ( ! (fp = fopen (outfile, "w")))
      fp = stderr;

    fprintf (fp, "$Id: gptl.c,v 1.153 2010/04/13 21:41:50 rosinski Exp $\n");
    fprintf (fp, "'count' is cumulative. All other stats are max/min\n");

    /* Print heading */

    fprintf (fp, "name");
    extraspace = max_name_len[0] - strlen ("name");
    for (n = 0; n < extraspace; ++n)
      fprintf (fp, " ");
    fprintf (fp, " count      wallmax (proc   thrd  )   wallmin (proc   thrd  )");

    for (n = 0; n < nevents; ++n) {
      fprintf (fp, " %8.8smax (proc   thrd  )", eventlist[n].str8);
      fprintf (fp, " %8.8smin (proc   thrd  )", eventlist[n].str8);
    }

    fprintf (fp, "\n");

    /*
    ** Gather and print stats based on list of thread 0 timers
    */

    for (ptr = timers[0]->next; ptr; ptr = ptr->next) {

      /* First, master gathers his own stats */

      get_threadstats (ptr->name, &summarystats);

      /* Broadcast a message to slaves asking for their results for this timer */

      ret = MPI_Bcast (ptr->name, MAX_CHARS+1, MPI_CHAR, 0, comm);

      /* Loop over slaves, receiving and processing the results */

      for (p = 1; p < nproc; ++p) {
	ret = MPI_Recv (&summarystats_slave, nbytes, MPI_BYTE, p, tag, comm, &status);
	if (summarystats_slave.count > 0)   /* timer found in slave */
	  get_summarystats (&summarystats, &summarystats_slave, p);
      }

      /* Print the results for this timer */

      fprintf (fp, "%s", ptr->name);
      extraspace = max_name_len[0] - strlen (ptr->name);
      for (n = 0; n < extraspace; ++n)
	fprintf (fp, " ");
      fprintf (fp, " %8lu %9.3f (%6d %6d) %9.3f (%6d %6d)", 
	       summarystats.count, 
	       summarystats.wallmax, summarystats.wallmax_p, summarystats.wallmax_t, 
	       summarystats.wallmin, summarystats.wallmin_p, summarystats.wallmin_t);
#ifdef HAVE_PAPI
      for (n = 0; n < nevents; ++n) {
	fprintf (fp, " %8.2e    (%6d %6d)", 
		 summarystats.papimax[n], summarystats.papimax_p[n], 
		 summarystats.papimax_t[n]);

	fprintf (fp, " %8.2e    (%6d %6d)", 
		 summarystats.papimin[n], summarystats.papimin_p[n], 
		 summarystats.papimin_t[n]);
      }
#endif
      fprintf (fp, "\n");
    }

    /* Signal that we're done */

    ret = MPI_Bcast (emptystring, MAX_CHARS+1, MPI_CHAR, 0, comm);

  } else {   /* iam != 0 (slave) */

    /* Loop until null message received from master */

    while (1) {
      ret = MPI_Bcast (name, MAX_CHARS+1, MPI_CHAR, 0, comm);
      if (strcmp (name, "") == 0)
	break;
      get_threadstats (name, &summarystats_slave);
      ret = MPI_Send (&summarystats_slave, nbytes, MPI_BYTE, 0, tag, comm);
    }
  }
  if (iam == 0 && fclose (fp) != 0)
    fprintf (stderr, "Attempt to close %s failed\n", outfile);
  return 0;
}

/* 
** get_threadstats: gather stats for timer "name" over all threads
**
** Input arguments:
**   name:  timer name
** Output arguments:
**   summarystats: max/min stats over all threads
*/

void get_threadstats (const char *name, 
		      Summarystats *summarystats)
{
#ifdef HAVE_PAPI
  int n;       /* event index */
#endif
  int t;       /* thread index */
  int indx;    /* returned from getentry() */
  Timer *ptr;  /* timer */

  /*
  ** This memset fortuitiously initializes the process values (_p) to master (0)
  */

  memset (summarystats, 0, sizeof (Summarystats));

  for (t = 0; t < GPTLnthreads; ++t) {
    if ((ptr = getentry (hashtable[t], name, &indx))) {

      summarystats->count += ptr->count;

      if (ptr->wall.accum > summarystats->wallmax) {
	summarystats->wallmax   = ptr->wall.accum;
	summarystats->wallmax_t = t;
      }

      if (ptr->wall.accum < summarystats->wallmin || summarystats->wallmin == 0.) {
	summarystats->wallmin   = ptr->wall.accum;
	summarystats->wallmin_t = t;
      }
#ifdef HAVE_PAPI
      for (n = 0; n < nevents; ++n) {
	double value;
	if (GPTL_PAPIget_eventvalue (eventlist[n].namestr, &ptr->aux, &value) != 0) {
	  fprintf (stderr, "Bad return from GPTL_PAPIget_eventvalue\n");
	  return;
	}
	if (value > summarystats->papimax[n]) {
	  summarystats->papimax[n]   = value;
	  summarystats->papimax_t[n] = t;
	}
	
	if (value < summarystats->papimin[n] || summarystats->papimin[n] == 0.) {
	  summarystats->papimin[n]   = value;
	  summarystats->papimin_t[n] = t;
	}
      }
#endif
    }
  }
}

/* 
** get_summarystats: write max/min stats into mpistats based on comparison
**                   with  summarystats_slave
**
** Input arguments:
**   p:              rank of summarystats_slave
**   summarystats_slave: stats from a slave process
** Input/Output arguments:
**   summarystats:       stats (starts out as master stats)
*/

void get_summarystats (Summarystats *summarystats, 
		       const Summarystats *summarystats_slave,
		       const int p)
{
  summarystats->count += summarystats_slave->count;

  if (summarystats_slave->wallmax > summarystats->wallmax) {
    summarystats->wallmax   = summarystats_slave->wallmax;
    summarystats->wallmax_p = p;
    summarystats->wallmax_t = summarystats_slave->wallmax_t;
  }
    
  if (summarystats_slave->wallmin < summarystats->wallmin || summarystats->wallmin == 0.) {
    summarystats->wallmin   = summarystats_slave->wallmin;
    summarystats->wallmin_p = p;
    summarystats->wallmin_t = summarystats_slave->wallmin_t;
  }
  
#ifdef HAVE_PAPI
  {
    int n;
    for (n = 0; n < nevents; ++n) {
      if (summarystats_slave->papimax[n] > summarystats->papimax[n]) {
	summarystats->papimax[n]   = summarystats_slave->papimax[n];
	summarystats->papimax_p[n] = p;
	summarystats->papimax_t[n] = summarystats_slave->papimax_t[n];
      }

      if (summarystats_slave->papimin[n] < summarystats->papimin[n] || 
	  summarystats->papimin[n] == 0.) {
	summarystats->papimin[n]   = summarystats_slave->papimin[n];
	summarystats->papimin_p[n] = p;
	summarystats->papimin_t[n] = summarystats_slave->papimin_t[n];
      }
    }
  }
#endif
}

/* 
** GPTLbarrier: When MPI enabled, set and time an MPI barrier
**
** Input arguments:
**   comm: commuicator (e.g. MPI_COMM_WORLD). If zero, use MPI_COMM_WORLD
**   name: region name
**
** Return value: 0 (success)
*/

int GPTLbarrier (MPI_Comm comm, const char *name)
{
  int ret;

  ret = GPTLstart (name);
  ret = MPI_Barrier (comm);
  ret = GPTLstop (name);
  return 0;
}
#endif

/*
** get_cpustamp: Invoke the proper system timer and return stats.
**
** Output arguments:
**   usr: user time
**   sys: system time
**
** Return value: 0 (success)
*/

static inline int get_cpustamp (long *usr, long *sys)
{
#ifdef HAVE_TIMES
  struct tms buf;

  (void) times (&buf);
  *usr = buf.tms_utime;
  *sys = buf.tms_stime;
  return 0;
#else
  return GPTLerror ("get_cpustamp: times() not available\n");
#endif
}

/*
** GPTLquery: return current status info about a timer. If certain stats are not 
** enabled, they should just have zeros in them. If PAPI is not enabled, input
** counter info is ignored.
** 
** Input args:
**   name:        timer name
**   maxcounters: max number of PAPI counters to get info for
**   t:           thread number (if < 0, the request is for the current thread)
**
** Output args:
**   count:            number of times this timer was called
**   onflg:            whether timer is currently on
**   wallclock:        accumulated wallclock time
**   usr:              accumulated user CPU time
**   sys:              accumulated system CPU time
**   papicounters_out: accumulated PAPI counters
*/

int GPTLquery (const char *name, 
	       int t,
	       int *count,
	       int *onflg,
	       double *wallclock,
	       double *dusr,
	       double *dsys,
	       long long *papicounters_out,
	       const int maxcounters)
{
  Timer *ptr;                /* linked list pointer */
  int indx;
  
  if ( ! initialized)
    return GPTLerror ("GPTLquery: GPTLinitialize has not been called\n");
  
  /*
  ** If t is < 0, assume the request is for the current thread
  */
  
  if (t < 0) {
    if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
      return GPTLerror ("GPTLquery: get_thread_num failure\n");
  } else {
    if (t >= maxthreads)
      return GPTLerror ("GPTLquery: requested thread %d is too big\n", t);
  }
  
  ptr = getentry (hashtable[t], name, &indx);
  if ( !ptr)
    return GPTLerror ("GPTLquery: requested timer %s does not exist\n", name);

  *onflg     = ptr->onflg;
  *count     = ptr->count;
  *wallclock = ptr->wall.accum;
  *dusr      = ptr->cpu.accum_utime / (double) ticks_per_sec;
  *dsys      = ptr->cpu.accum_stime / (double) ticks_per_sec;
#ifdef HAVE_PAPI
  GPTL_PAPIquery (&ptr->aux, papicounters_out, maxcounters);
#endif
  return 0;
}

/*
** GPTLquerycounters: return current PAPI counters for a timer.
** THIS ROUTINE ID DEPRECATED. USE GPTLget_eventvalue() instead
** 
** Input args:
**   name: timer name
**   t:    thread number (if < 0, the request is for the current thread)
**
** Output args:
**   papicounters_out: accumulated PAPI counters
*/

int GPTLquerycounters (const char *name, 
		       int t,
		       long long *papicounters_out)
{
  Timer *ptr;   /* linked list pointer */
  int indx;     /* hash index returned from getentry */
  
  if ( ! initialized)
    return GPTLerror ("GPTLquerycounters: GPTLinitialize has not been called\n");
  
  /*
  ** If t is < 0, assume the request is for the current thread
  */
  
  if (t < 0) {
    if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
      return GPTLerror ("GPTLquerycounters: get_thread_num failure\n");
  } else {
    if (t >= maxthreads)
      return GPTLerror ("GPTLquerycounters: requested thread %d is too big\n", t);
  }
  
  ptr = getentry (hashtable[t], name, &indx);
  if ( !ptr)
    return GPTLerror ("GPTLquerycounters: requested timer %s does not exist\n", name);

#ifdef HAVE_PAPI
  /* The 999 is a hack to say "give me all the counters" */
  GPTL_PAPIquery (&ptr->aux, papicounters_out, 999);
#endif
  return 0;
}

/*
** GPTLget_wallclock: return wallclock accumulation for a timer.
** 
** Input args:
**   timername: timer name
**   t:         thread number (if < 0, the request is for the current thread)
**
** Output args:
**   value: current wallclock accumulation for the timer
*/

int GPTLget_wallclock (const char *timername,
		      int t,
		      double *value)
{
  Timer *ptr; /* linked list pointer */
  int indx;   /* hash index returned from getentry (unused) */
  
  if ( ! initialized)
    return GPTLerror ("GPTLquery_event: GPTLinitialize has not been called\n");

  if ( ! wallstats.enabled)
    return GPTLerror ("GPTLquery_event: wallstats not enabled\n");
  
  /*
  ** If t is < 0, assume the request is for the current thread
  */
  
  if (t < 0) {
    if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
      return GPTLerror ("GPTLget_wallclock: get_thread_num failure\n");
  } else {
    if (t >= maxthreads)
      return GPTLerror ("GPTLget_wallclock: requested thread %d is too big\n", t);
  }
  
  ptr = getentry (hashtable[t], timername, &indx);
  if ( !ptr)
    return GPTLerror ("GPTLget_wallclock: requested timer %s does not exist\n", timername);

  *value = ptr->wall.accum;
  return 0;
}

/*
** GPTLget_eventvalue: return PAPI-based event value for a timer. All values will be
**   returned as doubles, even if the event is not derived.
** 
** Input args:
**   timername: timer name
**   eventname: event name (must be currently enabled)
**   t:         thread number (if < 0, the request is for the current thread)
**
** Output args:
**   value: current value of the event for this timer
*/

int GPTLget_eventvalue (const char *timername,
			const char *eventname,
			int t,
			double *value)
{
  Timer *ptr; /* linked list pointer */
  int indx;   /* hash index returned from getentry (unused) */
  
  if ( ! initialized)
    return GPTLerror ("GPTLquery_event: GPTLinitialize has not been called\n");
  
  /*
  ** If t is < 0, assume the request is for the current thread
  */
  
  if (t < 0) {
    if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
      return GPTLerror ("GPTLget_eventvalue: get_thread_num failure\n");
  } else {
    if (t >= maxthreads)
      return GPTLerror ("GPTLget_eventvalue: requested thread %d is too big\n", t);
  }
  
  ptr = getentry (hashtable[t], timername, &indx);
  if ( !ptr)
    return GPTLerror ("GPTLget_eventvalue: requested timer %s does not exist\n", timername);

#ifdef HAVE_PAPI
  return GPTL_PAPIget_eventvalue (eventname, &ptr->aux, value);
#else
  return GPTLerror ("GPTLquery_event: PAPI not enabled\n"); 
#endif
}

/*
** GPTLget_nregions: return number of regions (i.e. timer names) for this thread
** 
** Input args:
**   t:    thread number (if < 0, the request is for the current thread)
**
** Output args:
**   nregions: number of regions
*/

int GPTLget_nregions (int t, 
		      int *nregions)
{
  Timer *ptr;

  if ( ! initialized)
    return GPTLerror ("GPTLget_nregions: GPTLinitialize has not been called\n");
  
  /*
  ** If t is < 0, assume the request is for the current thread
  */
  
  if (t < 0) {
    if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
      return GPTLerror ("GPTLget_nregions: get_thread_num failure\n");
  } else {
    if (t >= maxthreads)
      return GPTLerror ("GPTLget_nregions: requested thread %d is too big\n", t);
  }
  
  *nregions = 0;
  for (ptr = timers[t]->next; ptr; ptr = ptr->next) 
    ++*nregions;

  return 0;
}

/*
** GPTLget_regionname: return region name for this thread
** 
** Input args:
**   t:      thread number (if < 0, the request is for the current thread)
**   region: region number
**   nc:     max number of chars to put in name
**
** Output args:
**   name    region name
*/

int GPTLget_regionname (int t,      /* thread number */
			int region, /* region number (0-based) */
			char *name, /* output region name */
			int nc)     /* number of chars in name (free from Fortran) */
{
  int ncpy;    /* number of characters to copy */
  int i;       /* index */
  Timer *ptr;

  if ( ! initialized)
    return GPTLerror ("GPTLget_nregionname: GPTLinitialize has not been called\n");
  
  /*
  ** If t is < 0, assume the request is for the current thread
  */
  
  if (t < 0) {
    if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0)
      return GPTLerror ("GPTLget_regionname: get_thread_num failure\n");
  } else {
    if (t >= maxthreads)
      return GPTLerror ("GPTLget_regionname: requested thread %d is too big\n", t);
  }
  
  ptr = timers[t]->next;
  for (i = 0; i < region; i++) {
    if ( ! ptr)
      return GPTLerror ("GPTLget_regionname: timer %d does not exist in thread %d\n",
			t, region);
    ptr = ptr->next;
  }

  if (ptr) {
    ncpy = MIN (nc, strlen (ptr->name));
    strncpy (name, ptr->name, ncpy);
    
    /*
    ** Adding the \0 is only important when called from C
    */

    if (ncpy < nc)
      name[ncpy] = '\0';
  } else {
    return GPTLerror ("GPTLget_regionname: timer %d does not exist in thread %d\n",
		      t, region);
  }
  return 0;
}

/*
** GPTLis_initialized: Return whether GPTL has been initialized
*/

int GPTLis_initialized (void)
{
  return (int) initialized;
}

/*
** GPTLpr_file_has_been_called: Return whether GPTLpr_file has been called
*/

int GPTLpr_has_been_called (void)
{
  return (int) pr_has_been_called;
}

/*
** getentry_instr: find hash table entry and return a pointer to it
**
** Input args:
**   hashtable: the hashtable (array)
**   self:      input address (from -finstrument-functions)
** Output args:
**   indx:      hashtable index
**
** Return value: pointer to the entry, or NULL if not found
*/

static inline Timer *getentry_instr (const Hashentry *hashtable, /* hash table */
				     void *self,                 /* address */
				     int *indx)                  /* hash index */
{
  int i;
  Timer *ptr = 0;

  /*
  ** Hash index is timer address modulo the table size
  ** On most machines, right-shifting the address helps because linkers often
  ** align functions on even boundaries
  */

  *indx = (((unsigned long) self) >> 4) % tablesize;
  for (i = 0; i < hashtable[*indx].nument; ++i) {
    if (hashtable[*indx].entries[i]->address == self) {
      ptr = hashtable[*indx].entries[i];
      break;
    }
  }
  return ptr;
}

/*
** getentry: find the entry in the hash table and return a pointer to it.
**
** Input args:
**   hashtable: the hashtable (array)
**   name:      string to be hashed on (specifically, summed)
** Output args:
**   indx:      hashtable index
**
** Return value: pointer to the entry, or NULL if not found
*/

static inline Timer *getentry (const Hashentry *hashtable, /* hash table */
			       const char *name,           /* name to hash */
			       int *indx)                  /* hash index */
{
  int i;                      /* loop index */
  const unsigned char *c;     /* pointer to elements of "name" */
  Timer *ptr = 0;

  /* Hash value is sum of chars times their position index */

  *indx = 0;
  c = (unsigned char *) name;
  for (i = 1; *c && i < MAX_CHARS+1; ++c, ++i) {
    *indx += (*c) * i;
  }

  *indx %= tablesize;

  /* 
  ** If nument exceeds 1 there was a hash collision and we must search
  ** linearly through an array for a match
  */

  for (i = 0; i < hashtable[*indx].nument; i++) {
    if (STRMATCH (name, hashtable[*indx].entries[i]->name)) {
      ptr = hashtable[*indx].entries[i];
      break;
    }
  }
  return ptr;
}

/*
** Add entry points for auto-instrumented codes
** (-finstrument-functions in gcc, pathcc,
**  -Minstrument:functions in pgcc)
*/

#ifdef __cplusplus
extern "C" {
#endif

void __cyg_profile_func_enter (void *this_fn,
                               void *call_site)
{
  (void) GPTLstart_instr (this_fn);
}

void __cyg_profile_func_exit (void *this_fn,
                              void *call_site)
{
  (void) GPTLstop_instr (this_fn);
}

#ifdef __cplusplus
};
#endif

#ifdef HAVE_NANOTIME
#ifdef BIT64
/* 64-bit code copied from PAPI library */
static inline unsigned long long nanotime (void)
{
    unsigned long long val;
    do {
      unsigned int a,d;
      asm volatile("rdtsc" : "=a" (a), "=d" (d));
      (val) = ((unsigned long)a) | (((unsigned long)d)<<32);
    } while(0);

    return (val);
}
#else
static inline unsigned long long nanotime (void)
{
  unsigned long long val;
  __asm__ __volatile__("rdtsc" : "=A" (val) : );
  return (val);
}
#endif

#define LEN 4096

static float get_clockfreq ()
{
  FILE *fd = 0;
  char buf[LEN];
  int is;

  if ( ! (fd = fopen ("/proc/cpuinfo", "r"))) {
    fprintf (stderr, "get_clockfreq: can't open /proc/cpuinfo\n");
    return -1.;
  }

  while (fgets (buf, LEN, fd)) {
    if (strncmp (buf, "cpu MHz", 7) == 0) {
      for (is = 7; buf[is] != '\0' && !isdigit (buf[is]); is++);
      if (isdigit (buf[is]))
	return atof (&buf[is]);
    }
  }

  return -1.;
}
#endif

/*
** The following are the set of underlying timing routines which may or may
** not be available. And their accompanying init routines.
** NANOTIME is currently only available on x86.
*/

static int init_nanotime ()
{
#ifdef HAVE_NANOTIME
  if ((cpumhz = get_clockfreq ()) < 0)
    return GPTLerror ("Can't get clock freq\n");

  if (verbose)
    printf ("init_nanotime: Clock rate = %f MHz\n", cpumhz);
  cyc2sec = 1./(cpumhz * 1.e6);
  return 0;
#else
  return GPTLerror ("init_nanotime: not enabled\n");
#endif
}

static inline double utr_nanotime ()
{
#ifdef HAVE_NANOTIME
  double timestamp;
  timestamp = nanotime () * cyc2sec;
  return timestamp;
#else
  (void) GPTLerror ("utr_nanotime: not enabled\n");
  return -1.;
#endif
}

/*
** rtc is currently only available on UNICOSMP
*/

static int init_rtc ()
{
#ifdef UNICOSMP
  extern long long rtc_rate_();
  ticks2sec = 1./rtc_rate_();
  if (verbose)
    printf ("init_rtc: ticks per sec=%g\n", rtc_rate_());
  return 0;
#else
  return GPTLerror ("init_rtc: not enabled\n");
#endif
}
  
static inline double utr_rtc ()
{
#ifdef UNICOSMP
  return _rtc () * ticks2sec;
#else
  (void) GPTLerror ("utr_rtc: not enabled\n");
  return -1.;
#endif
}

/*
** MPI_Wtime requires the MPI lib.
*/

static int init_mpiwtime ()
{
#ifdef HAVE_MPI
  return 0;
#else
  return GPTLerror ("utr_mpiwtime: not enabled\n");
#endif
}

static inline double utr_mpiwtime ()
{
#ifdef HAVE_MPI
  return MPI_Wtime ();
#else
  (void) GPTLerror ("utr_mpiwtime: not enabled\n");
  return -1.;
#endif
}

/*
** PAPI_get_real_usec requires the PAPI lib.
*/

static int init_papitime ()
{
#ifdef HAVE_PAPI
  ref_papitime = PAPI_get_real_usec ();
  if (verbose)
    printf ("init_papitime: ref_papitime=%ld\n", (long) ref_papitime);
  return 0;
#else
  return GPTLerror ("init_papitime: not enabled\n");
#endif
}
  
static inline double utr_papitime ()
{
#ifdef HAVE_PAPI
  return (PAPI_get_real_usec () - ref_papitime) * 1.e-6;
#else
  (void) GPTLerror ("utr_papitime: not enabled\n");
  return -1.;
#endif
}

/* 
** Probably need to link with -lrt for this one to work 
*/

static int init_clock_gettime ()
{
#ifdef HAVE_LIBRT
  struct timespec tp;
  (void) clock_gettime (CLOCK_REALTIME, &tp);
  ref_clock_gettime = tp.tv_sec;
  if (verbose)
    printf ("init_clock_gettime: ref_clock_gettime=%ld\n", (long) ref_clock_gettime);
  return 0;
#else
  return GPTLerror ("init_clock_gettime: not enabled\n");
#endif
}

static inline double utr_clock_gettime ()
{
#ifdef HAVE_LIBRT
  struct timespec tp;
  (void) clock_gettime (CLOCK_REALTIME, &tp);
  return (tp.tv_sec - ref_clock_gettime) + 1.e-9*tp.tv_nsec;
#else
  (void) GPTLerror ("utr_clock_gettime: not enabled\n");
  return -1.;
#endif
}

/*
** Default available most places: gettimeofday
*/

static int init_gettimeofday ()
{
#ifdef HAVE_GETTIMEOFDAY
  struct timeval tp;
  (void) gettimeofday (&tp, 0);
  ref_gettimeofday = tp.tv_sec;
  if (verbose)
    printf ("init_gettimeofday: ref_gettimeofday=%ld\n", (long) ref_gettimeofday);
  return 0;
#else
  return GPTLerror ("init_gettimeofday: not enabled\n");
#endif
}

static inline double utr_gettimeofday ()
{
#ifdef HAVE_GETTIMEOFDAY
  struct timeval tp;
  (void) gettimeofday (&tp, 0);
  return (tp.tv_sec - ref_gettimeofday) + 1.e-6*tp.tv_usec;
#else
  return GPTLerror ("utr_gettimeofday: not enabled\n");
#endif
}

/* 
** Determine underlying timing routine overhead: call it 100 times.
*/

static double utr_getoverhead ()
{
  double val1;
  double val2;
  int i;

  val1 = (*ptr2wtimefunc)();
  for (i = 0; i < 10; ++i) {
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
    val2 = (*ptr2wtimefunc)();
  }
  return 0.01 * (val2 - val1);
}

/*
** printself_andchildren: Recurse through call tree, printing stats for self, then children
*/

static void printself_andchildren (const Timer *ptr,
				   FILE *fp, 
				   const int t,
				   const int depth,
				   const double tot_overhead)
{
  int n;

  if (depth > -1)     /* -1 flag is to avoid printing stats for dummy outer timer */
    printstats (ptr, fp, t, depth, true, tot_overhead);

  for (n = 0; n < ptr->nchildren; n++)
    printself_andchildren (ptr->children[n], fp, t, depth+1, tot_overhead);
}

#ifdef ENABLE_PMPI
Timer *GPTLgetentry (char *name)
{
  int t;
  int indx;
  
  if ( ! initialized) {
    (void) GPTLerror ("GPTLgetentry: initialization was not completed\n");
    return 0;
  }

  if ((t = get_thread_num (&GPTLnthreads, &maxthreads)) < 0) {
    (void) GPTLerror ("GPTLgetentry: bad return from get_thread_num\n");
    return 0;
  }

  return (getentry (hashtable[t], name, &indx));
}
#endif
