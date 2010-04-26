/*
** $Id: gptl.h,v 1.56 2010/02/18 04:23:08 rosinski Exp $
**
** Author: Jim Rosinski
**
** #include this file in C and C++ codes which invoke GPTL functions.
*/

#ifndef GPTL_H
#define GPTL_H

/*
** Options settable by a call to GPTLsetoption() (default in parens)
** These numbers need to be small integers because GPTLsetoption can
** be passed PAPI counters, and we need to avoid collisions in that
** integer space. PAPI presets are big negative integers, and PAPI
** native events are big positive integers.
*/

typedef enum {
  GPTLsync_mpi       = 0,  /* Synchronize before certain MPI calls (PMPI-mode only) */
  GPTLwall           = 1,  /* Collect wallclock stats (true) */
  GPTLcpu            = 2,  /* Collect CPU stats (false)*/
  GPTLabort_on_error = 3,  /* Abort on failure (false) */
  GPTLoverhead       = 4,  /* Estimate overhead of underlying timing routine (true) */
  GPTLdepthlimit     = 5,  /* Only print timers this depth or less in the tree (inf) */
  GPTLverbose        = 6,  /* Verbose output (false) */
  GPTLnarrowprint    = 7,  /* Print PAPI and derived stats in 8 columns not 16 (true) */
  GPTLpercent        = 9,  /* Add a column for percent of first timer (false) */
  GPTLpersec         = 10, /* Add a PAPI column that prints "per second" stats (true) */
  GPTLmultiplex      = 11, /* Allow PAPI multiplexing (false) */
  GPTLdopr_preamble  = 12, /* Print preamble info (true) */
  GPTLdopr_threadsort= 13, /* Print sorted thread stats (true) */
  GPTLdopr_multparent= 14, /* Print multiple parent info (true) */
  GPTLdopr_collision = 15, /* Print hastable collision info (true) */
  GPTLprint_method   = 16, /* Tree print method: first parent, last parent
			      most frequent, or full tree (most frequent) */
  /*
  ** These are derived counters based on PAPI counters. All default to false
  */
  GPTL_IPC           = 17, /* Instructions per cycle */
  GPTL_CI            = 18, /* Computational intensity */
  GPTL_FPC           = 19, /* FP ops per cycle */
  GPTL_FPI           = 20, /* FP ops per instruction */
  GPTL_LSTPI         = 21, /* Load-store instruction fraction */
  GPTL_DCMRT         = 22, /* L1 miss rate (fraction) */
  GPTL_LSTPDCM       = 23, /* Load-stores per L1 miss */
  GPTL_L2MRT         = 24, /* L2 miss rate (fraction) */
  GPTL_LSTPL2M       = 25, /* Load-stores per L2 miss */
  GPTL_L3MRT         = 26  /* L3 read miss rate (fraction) */
} Option;

/*
** Underlying wallclock timer: optimize for best granularity with least overhead.
** These numbers need not be distinct from the above because these are passed
** to GPTLsetutr() and the above are passed to GPTLsetoption()
*/

typedef enum {
  GPTLgettimeofday   = 1, /* the default */
  GPTLnanotime       = 2, /* only available on x86 */
  GPTLrtc            = 3, /* real time clock (UNICOSMP only) */
  GPTLmpiwtime       = 4, /* MPI_Wtime */
  GPTLclockgettime   = 5, /* clock_gettime */
  GPTLpapitime       = 6  /* only if PAPI is available */
} Funcoption;

/*
** How to report parent/child relationships at print time (for children with multiple parents)
*/

typedef enum {
  GPTLfirst_parent  = 1,  /* first parent found */
  GPTLlast_parent   = 2,  /* last parent found */
  GPTLmost_frequent = 3,  /* most frequent parent (default) */
  GPTLfull_tree     = 4   /* complete call tree */
} Method;

/*
** Function prototypes
*/

#ifdef __cplusplus
extern "C" {
#endif

extern int GPTLsetoption (const int, const int);
extern int GPTLinitialize (void);
extern int GPTLstart (const char *);
extern int GPTLstop (const char *);
extern int GPTLstamp (double *, double *, double *);
extern int GPTLpr (const int);
extern int GPTLpr_file (const char *);

/*
** Use K&R prototype for these 2 because they require MPI
*/
extern int GPTLpr_summary ();
extern int GPTLbarrier ();

extern int GPTLreset (void);
extern int GPTLfinalize (void);
extern int GPTLget_memusage (int *, int *, int *, int *, int *);
extern int GPTLprint_memusage (const char *);
extern int GPTLenable (void);
extern int GPTLdisable (void);
extern int GPTLsetutr (const int);
extern int GPTLquery (const char *, int, int *, int *, double *, double *, double *,
		      long long *, const int);
extern int GPTLquerycounters (const char *, int, long long *);
extern int GPTLget_wallclock (const char *, int, double *);
extern int GPTLget_eventvalue (const char *, const char *, int, double *);
extern int GPTLget_nregions (int, int *);
extern int GPTLget_regionname (int, int, char *, int);
extern int GPTL_PAPIlibraryinit (void);
extern int GPTLevent_name_to_code (const char *, int *);
extern int GPTLevent_code_to_name (const int, char *);

#ifdef __cplusplus
};
#endif

#endif
