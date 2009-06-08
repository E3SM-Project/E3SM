#ifdef HAVE_PAPI

#include <papi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "private.h"
#include "gptl.h"

#if ( defined THREADED_OMP )
#include <omp.h>
#elif ( defined THREADED_PTHREADS )
#include <pthread.h>
#endif

/* Mapping of PAPI counters to short and long printed strings */

static const Entry papitable [] = {
  {PAPI_L1_DCM, "PAPI_L1_DCM", "L1 Dcache miss  ", "Level 1 data cache misses"},
  {PAPI_L1_ICM, "PAPI_L1_ICM", "L1 Icache miss  ", "Level 1 instruction cache misses"},
  {PAPI_L2_DCM, "PAPI_L2_DCM", "L2 Dcache miss  ", "Level 2 data cache misses"},
  {PAPI_L2_ICM, "PAPI_L2_ICM", "L2 Icache miss  ", "Level 2 instruction cache misses"},
  {PAPI_L3_DCM, "PAPI_L3_DCM", "L3 Dcache miss  ", "Level 3 data cache misses"},
  {PAPI_L3_ICM, "PAPI_L3_ICM", "L3 Icache miss  ", "Level 3 instruction cache misses"},
  {PAPI_L1_TCM, "PAPI_L1_TCM", "L1 cache miss   ", "Level 1 total cache misses"},
  {PAPI_L2_TCM, "PAPI_L2_TCM", "L2 cache miss   ", "Level 2 total cache misses"},
  {PAPI_L3_TCM, "PAPI_L3_TCM", "L3 cache miss   ", "Level 3 total cache misses"},
  {PAPI_CA_SNP, "PAPI_CA_SNP", "Snoops          ", "Snoops          "},
  {PAPI_CA_SHR, "PAPI_CA_SHR", "PAPI_CA_SHR     ", "Request for shared cache line (SMP)"},
  {PAPI_CA_CLN, "PAPI_CA_CLN", "PAPI_CA_CLN     ", "Request for clean cache line (SMP)"},
  {PAPI_CA_INV, "PAPI_CA_INV", "PAPI_CA_INV     ", "Request for cache line Invalidation (SMP)"},
  {PAPI_CA_ITV, "PAPI_CA_ITV", "PAPI_CA_ITV     ", "Request for cache line Intervention (SMP)"},
  {PAPI_L3_LDM, "PAPI_L3_LDM", "L3 load misses  ", "Level 3 load misses"},
  {PAPI_L3_STM, "PAPI_L3_STM", "L3 store misses ", "Level 3 store misses"},
  {PAPI_BRU_IDL,"PAPI_BRU_IDL","PAPI_BRU_IDL    ", "Cycles branch units are idle"},
  {PAPI_FXU_IDL,"PAPI_FXU_IDL","PAPI_FXU_IDL    ", "Cycles integer units are idle"},
  {PAPI_FPU_IDL,"PAPI_FPU_IDL","PAPI_FPU_IDL    ", "Cycles floating point units are idle"},
  {PAPI_LSU_IDL,"PAPI_LSU_IDL","PAPI_LSU_IDL    ", "Cycles load/store units are idle"},
  {PAPI_TLB_DM,"PAPI_TLB_DM"  ,"Data TLB misses ", "Data translation lookaside buffer misses"},
  {PAPI_TLB_IM,"PAPI_TLB_IM", "Inst TLB misses ", "Instr translation lookaside buffer misses"},
  {PAPI_TLB_TL,"PAPI_TLB_TL", "Tot TLB misses  ", "Total translation lookaside buffer misses"},
  {PAPI_L1_LDM,"PAPI_L1_LDM", "L1 load misses  ", "Level 1 load misses"},
  {PAPI_L1_STM,"PAPI_L1_STM", "L1 store misses ", "Level 1 store misses"},
  {PAPI_L2_LDM,"PAPI_L2_LDM", "L2 load misses  ", "Level 2 load misses"},
  {PAPI_L2_STM,"PAPI_L2_STM", "L2 store misses ", "Level 2 store misses"},
  {PAPI_BTAC_M,"PAPI_BTAC_M", "BTAC miss       ", "BTAC miss"},
  {PAPI_PRF_DM,"PAPI_PRF_DM", "PAPI_PRF_DM     ", "Prefetch data instruction caused a miss"},
  {PAPI_L3_DCH,"PAPI_L3_DCH", "L3 DCache Hit   ", "Level 3 Data Cache Hit"},
  {PAPI_TLB_SD,"PAPI_TLB_SD", "PAPI_TLB_SD     ", "Xlation lookaside buffer shootdowns (SMP)"},
  {PAPI_CSR_FAL,"PAPI_CSR_FAL","PAPI_CSR_FAL    ", "Failed store conditional instructions"},
  {PAPI_CSR_SUC,"PAPI_CSR_SUC","PAPI_CSR_SUC    ", "Successful store conditional instructions"},
  {PAPI_CSR_TOT,"PAPI_CSR_TOT","PAPI_CSR_TOT    ", "Total store conditional instructions"},
  {PAPI_MEM_SCY,"PAPI_MEM_SCY","Cyc Stalled Mem ", "Cycles Stalled Waiting for Memory Access"},
  {PAPI_MEM_RCY,"PAPI_MEM_RCY","Cyc Stalled MemR", "Cycles Stalled Waiting for Memory Read"},
  {PAPI_MEM_WCY,"PAPI_MEM_WCY","Cyc Stalled MemW", "Cycles Stalled Waiting for Memory Write"},
  {PAPI_STL_ICY,"PAPI_STL_ICY","Cyc no InstrIss ", "Cycles with No Instruction Issue"},
  {PAPI_FUL_ICY,"PAPI_FUL_ICY","Cyc Max InstrIss", "Cycles with Maximum Instruction Issue"},
  {PAPI_STL_CCY,"PAPI_STL_CCY","Cyc No InstrComp", "Cycles with No Instruction Completion"},
  {PAPI_FUL_CCY,"PAPI_FUL_CCY","Cyc Max InstComp", "Cycles with Maximum Instruction Completion"},
  {PAPI_HW_INT,"PAPI_HW_INT", "HW interrupts   ", "Hardware interrupts"},
  {PAPI_BR_UCN,"PAPI_BR_UCN", "Uncond br instr ", "Unconditional branch instructions executed"},
  {PAPI_BR_CN,"PAPI_BR_CN",  "Cond br instr ex", "Conditional branch instructions executed"},
  {PAPI_BR_TKN,"PAPI_BR_TKN", "Cond br instr tk", "Conditional branch instructions taken"},
  {PAPI_BR_NTK,"PAPI_BR_NTK", "Cond br instrNtk", "Conditional branch instructions not taken"},
  {PAPI_BR_MSP,"PAPI_BR_MSP", "Cond br instrMPR", "Conditional branch instructions mispred"},
  {PAPI_BR_PRC,"PAPI_BR_PRC", "Cond br instrCPR", "Conditional branch instructions corr. pred"},
  {PAPI_FMA_INS,"PAPI_FMA_INS","FMA instr comp  ", "FMA instructions completed"},
  {PAPI_TOT_IIS,"PAPI_TOT_IIS","Total instr iss ", "Total instructions issued"},
  {PAPI_TOT_INS,"PAPI_TOT_INS","Total instr ex  ", "Total instructions executed"},
  {PAPI_INT_INS,"PAPI_INT_INS","Int instr ex    ", "Integer instructions executed"},
  {PAPI_FP_INS, "PAPI_FP_INS", "FP instr ex     ", "Floating point instructions executed"},
  {PAPI_LD_INS,"PAPI_LD_INS", "Load instr ex   ", "Load instructions executed"},
  {PAPI_SR_INS,"PAPI_SR_INS", "Store instr ex  ", "Store instructions executed"},
  {PAPI_BR_INS,"PAPI_BR_INS", "br instr ex     ", "Total branch instructions executed"},
  {PAPI_VEC_INS,"PAPI_VEC_INS","Vec/SIMD instrEx", "Vector/SIMD instructions executed"},
  {PAPI_RES_STL,"PAPI_RES_STL","Cyc proc stalled", "Cycles processor is stalled on resource"},
  {PAPI_FP_STAL,"PAPI_FP_STAL","Cyc any FP stall", "Cycles any FP units are stalled"},
  {PAPI_TOT_CYC,"PAPI_TOT_CYC","Total cycles    ", "Total cycles"},
  {PAPI_LST_INS,"PAPI_LST_INS","Tot L/S inst ex ", "Total load/store inst. executed"},
  {PAPI_SYC_INS,"PAPI_SYC_INS","Sync. inst. ex  ", "Sync. inst. executed"},
  {PAPI_L1_DCH,"PAPI_L1_DCH", "L1 D Cache Hit  ", "L1 D Cache Hit"},
  {PAPI_L2_DCH,"PAPI_L2_DCH", "L2 D Cache Hit  ", "L2 D Cache Hit"},
  {PAPI_L1_DCA,"PAPI_L1_DCA", "L1 D Cache Acc  ", "L1 D Cache Access"},
  {PAPI_L2_DCA,"PAPI_L2_DCA", "L2 D Cache Acc  ", "L2 D Cache Access"},
  {PAPI_L3_DCA,"PAPI_L3_DCA", "L3 D Cache Acc  ", "L3 D Cache Access"},
  {PAPI_L1_DCR,"PAPI_L1_DCR", "L1 D Cache Read ", "L1 D Cache Read"},
  {PAPI_L2_DCR,"PAPI_L2_DCR", "L2 D Cache Read ", "L2 D Cache Read"},
  {PAPI_L3_DCR,"PAPI_L3_DCR", "L3 D Cache Read ", "L3 D Cache Read"},
  {PAPI_L1_DCW,"PAPI_L1_DCW", "L1 D Cache Write", "L1 D Cache Write"},
  {PAPI_L2_DCW,"PAPI_L2_DCW", "L2 D Cache Write", "L2 D Cache Write"},
  {PAPI_L3_DCW,"PAPI_L3_DCW", "L3 D Cache Write", "L3 D Cache Write"},
  {PAPI_L1_ICH,"PAPI_L1_ICH", "L1 I cache hits ", "L1 instruction cache hits"},
  {PAPI_L2_ICH,"PAPI_L2_ICH", "L2 I cache hits ", "L2 instruction cache hits"},
  {PAPI_L3_ICH,"PAPI_L3_ICH", "L3 I cache hits ", "L3 instruction cache hits"},
  {PAPI_L1_ICA,"PAPI_L1_ICA", "L1 I cache acc  ", "L1 instruction cache accesses"},
  {PAPI_L2_ICA,"PAPI_L2_ICA", "L2 I cache acc  ", "L2 instruction cache accesses"},
  {PAPI_L3_ICA,"PAPI_L3_ICA", "L3 I cache acc  ", "L3 instruction cache accesses"},
  {PAPI_L1_ICR,"PAPI_L1_ICR", "L1 I cache reads", "L1 instruction cache reads"},
  {PAPI_L2_ICR,"PAPI_L2_ICR", "L2 I cache reads", "L2 instruction cache reads"},
  {PAPI_L3_ICR,"PAPI_L3_ICR", "L3 I cache reads", "L3 instruction cache reads"},
  {PAPI_L1_ICW,"PAPI_L1_ICW", "L1 I cache write", "L1 instruction cache writes"},
  {PAPI_L2_ICW,"PAPI_L2_ICW", "L2 I cache write", "L2 instruction cache writes"},
  {PAPI_L3_ICW,"PAPI_L3_ICW", "L3 I cache write", "L3 instruction cache writes"},
  {PAPI_L1_TCH,"PAPI_L1_TCH", "L1 cache hits   ", "L1 total cache hits"},
  {PAPI_L2_TCH,"PAPI_L2_TCH", "L2 cache hits   ", "L2 total cache hits"},
  {PAPI_L3_TCH,"PAPI_L3_TCH", "L3 cache hits   ", "L3 total cache hits"},
  {PAPI_L1_TCA,"PAPI_L1_TCA", "L1 cache access ", "L1 total cache accesses"},
  {PAPI_L2_TCA,"PAPI_L2_TCA", "L2 cache access ", "L2 total cache accesses"},
  {PAPI_L3_TCA,"PAPI_L3_TCA", "L3 cache access ", "L3 total cache accesses"},
  {PAPI_L1_TCR,"PAPI_L1_TCR", "L1 cache reads  ", "L1 total cache reads"},
  {PAPI_L2_TCR,"PAPI_L2_TCR", "L2 cache reads  ", "L2 total cache reads"},
  {PAPI_L3_TCR,"PAPI_L3_TCR", "L3 cache reads  ", "L3 total cache reads"},
  {PAPI_L1_TCW,"PAPI_L1_TCW", "L1 cache writes ", "L1 total cache writes"},
  {PAPI_L2_TCW,"PAPI_L2_TCW", "L2 cache writes ", "L2 total cache writes"},
  {PAPI_L3_TCW,"PAPI_L3_TCW", "L3 cache writes ", "L3 total cache writes"},
  {PAPI_FML_INS,"PAPI_FML_INS","FM ins          ", "FM ins"},
  {PAPI_FAD_INS,"PAPI_FAD_INS","FA ins          ", "FA ins"},
  {PAPI_FDV_INS,"PAPI_FDV_INS","FD ins          ", "FD ins"},
  {PAPI_FSQ_INS,"PAPI_FSQ_INS","FSq ins         ", "FSq ins"},
  {PAPI_FNV_INS,"PAPI_FNV_INS","Finv ins        ", "Finv ins"},
  {PAPI_FP_OPS, "PAPI_FP_OPS", "FP ops executed ", "Floating point operations executed"}
};

static const int npapientries = sizeof (papitable) / sizeof (Entry);
static int papieventlist[MAX_AUX];        /* list of PAPI events to be counted */
static Pr_event pr_event[MAX_AUX];        /* list of events (PAPI or derived) */

/* Derived events */
static const Entry derivedtable [] = {
  {GPTL_IPC, "GPTL_IPC", "Instr per cycle ", "Instructions per cycle"},
  {GPTL_CI,  "GPTL_CI",  "Comp Intensity  ", "Computational intensity"}
};
static const int nderivedentries = sizeof (derivedtable) / sizeof (Entry);

static int npapievents = 0;              /* number of PAPI events: initialize to 0 */ 
static int nevents = 0;                  /* number of events: initialize to 0 */ 
static int *EventSet;                    /* list of events to be counted by PAPI */
static long_long **papicounters;         /* counters returned from PAPI */

static char papiname[PAPI_MAX_STR_LEN];  /* returned from PAPI_event_code_to_name */
static const int BADCOUNT = -999999;     /* Set counters to this when they are bad */
static bool is_multiplexed = false;      /* whether multiplexed (always start false)*/
static bool narrowprint = true;          /* only use 8 digits not 16 for counter prints */
static bool persec = true;               /* print PAPI stats per second */
static bool enable_multiplexing = true;  /* whether to try multiplexing */
static bool verbose = false;             /* output verbosity */

/* Function prototypes */

static int create_and_start_events (const int);
static int canenable (int);
static int canenable2 (int, int);
static int is_enabled (int);
static int enable (int);
static int getderivedidx (int);

/*
** GPTL_PAPIsetoption: enable or disable PAPI event defined by "counter". Called 
**   from GPTLsetoption.  Since all events are off by default, val=false degenerates
**   to a no-op.  Coded this way to be consistent with the rest of GPTL
**
** Input args: 
**   counter: PAPI counter
**   val:     true or false for enable or disable
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTL_PAPIsetoption (const int counter,  /* PAPI counter (or option) */
			const int val)      /* true or false for enable or disable */
{
  int n;       /* loop index */
  int ret;     /* return code */
  int numidx;  /* numerator index */
  int idx;

  /*
  ** First, check for option which is not an actual counter
  */

  if (counter == GPTLmultiplex) {
    enable_multiplexing = (bool) val;
    if (verbose)
      printf ("GPTL_PAPIsetoption: set GPTLmultiplex to %d\n", val);
    return 0;
  }

  if (counter == GPTLnarrowprint) {
    narrowprint = (bool) val;
    if (verbose)
      printf ("GPTL_PAPIsetoption: set GPTLnarrowprint to %d\n", val);
    return 0;
  }

  if (counter == GPTLpersec) {
    persec = (bool) val;
    return 0;
  }

  /* Initialize PAPI if it hasn't already been done */

  if (GPTL_PAPIlibraryinit () < 0)
    return GPTLerror ("GPTL_PAPIsetoption: PAPI library init error\n");

  /* Ensure max nevents won't be exceeded */

  if (nevents+1 > MAX_AUX)
    return GPTLerror ("GPTL_PAPIsetoption: %d is too many events\n", nevents+1);

  /* Check derived events */

  switch (counter) {
  case GPTL_IPC:
    if ( ! canenable2 (PAPI_TOT_INS, PAPI_TOT_CYC))
      return GPTLerror ("GPTL_PAPIsetoption: canenable2 return says GPTL_IPC unavailable\n");

    idx = getderivedidx (GPTL_IPC);
    pr_event[nevents].event    = derivedtable[idx];
    pr_event[nevents].numidx   = enable (PAPI_TOT_INS);
    pr_event[nevents].denomidx = enable (PAPI_TOT_CYC);
    ++nevents;
    return 0;
  case GPTL_CI:
    if ( ! canenable2 (PAPI_FP_OPS, PAPI_LST_INS))
      return GPTLerror ("GPTL_PAPIsetoption: canenable2 return says GPTL_CI unavailable\n");

    idx = getderivedidx (GPTL_CI);
    pr_event[nevents].event    = derivedtable[idx];
    pr_event[nevents].numidx   = enable (PAPI_FP_OPS);
    pr_event[nevents].denomidx = enable (PAPI_LST_INS);
    ++nevents;
    return 0;
  default:
    break;
  }

  /* Check PAPI presets */

  for (n = 0; n < npapientries; n++) {
    if (counter == papitable[n].counter) {
      if ((numidx = is_enabled (counter)) >= 0) {
	pr_event[nevents].event  = papitable[n];
	pr_event[nevents].numidx = numidx;
      } else if (canenable (counter)) {
	pr_event[nevents].event  = papitable[n];
	pr_event[nevents].numidx = enable (counter);
      } else {
	return GPTLerror ("GPTL_PAPIsetoption: Can't enable event \n", 
			  papitable[n].str);
      }
      if (verbose)
	printf ("GPTL_PAPIsetoption: will print event %s\n", 
		pr_event[nevents].event.str);
      ++nevents;
      return 0;
    }
  }

  /*
  ** Check native events last: If PAPI_event_code_to_name fails, give up
  */
  
  if ((ret = PAPI_event_code_to_name (counter, papiname)) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIsetoption: PAPI_strerror: %s\n", PAPI_strerror (ret));

  /*
  ** A table with predefined names of various lengths does not exist for
  ** native events. Just truncate papiname.
  */

  if ((numidx = is_enabled (counter)) >= 0) {
    pr_event[nevents].event.counter    = counter;

    pr_event[nevents].event.counterstr = GPTLallocate (12+1);
    strncpy (pr_event[nevents].event.counterstr, papiname, 12);
    pr_event[nevents].event.counterstr[12] = '\0';

    pr_event[nevents].event.prstr      = GPTLallocate (16+1);
    strncpy (pr_event[nevents].event.prstr, papiname, 16);
    pr_event[nevents].event.prstr[16] = '\0';

    pr_event[nevents].event.str        = GPTLallocate (PAPI_MAX_STR_LEN);
    strncpy (pr_event[nevents].event.str, papiname, PAPI_MAX_STR_LEN);

    pr_event[nevents].numidx           = numidx;
  } else if (canenable (counter)) {
    pr_event[nevents].event.counter    = counter;

    pr_event[nevents].event.counterstr = GPTLallocate (12+1);
    strncpy (pr_event[nevents].event.counterstr, papiname, 12);
    pr_event[nevents].event.counterstr[12] = '\0';

    pr_event[nevents].event.prstr      = GPTLallocate (16+1);
    strncpy (pr_event[nevents].event.prstr, papiname, 16);
    pr_event[nevents].event.prstr[16] = '\0';

    pr_event[nevents].event.str        = GPTLallocate (PAPI_MAX_STR_LEN);
    strncpy (pr_event[nevents].event.str, papiname, PAPI_MAX_STR_LEN);

    pr_event[nevents].numidx           = enable (counter);
  } else {
    return GPTLerror ("GPTL_PAPIsetoption: Can't enable event %s\n", papiname);
  }

  if (verbose)
    printf ("GPTL_PAPIsetoption: will print event %s\n", pr_event[nevents].event.str);

  ++nevents;
  return 0;
}

/*
** canenable: determine whether a PAPI counter can be enabled
**
** Input args: 
**   counter: PAPI counter
**
** Return value: 0 (success) or non-zero (failure)
*/
 
int canenable (int counter)
{
  if (npapievents+1 > MAX_AUX)
    return false;

  if (PAPI_query_event (counter) != PAPI_OK) {
    (void) PAPI_event_code_to_name (counter, papiname);
    fprintf (stderr, "canenable: event %s not available on this arch\n", papiname);
    return false;
  }

  return true;
}

/*
** canenable2: determine whether 2 PAPI counters can be enabled
**
** Input args: 
**   counter1: PAPI counter
**   counter2: PAPI counter
**
** Return value: 0 (success) or non-zero (failure)
*/
 
int canenable2 (int counter1, int counter2)
{
  if (npapievents+2 > MAX_AUX)
    return false;

  if (PAPI_query_event (counter1) != PAPI_OK) {
    (void) PAPI_event_code_to_name (counter1, papiname);
    fprintf (stderr, "canenable2: event %s not available on this arch\n", papiname);
    return false;
  }

  if (PAPI_query_event (counter2) != PAPI_OK) {
    (void) PAPI_event_code_to_name (counter2, papiname);
    fprintf (stderr, "canenable2: event %s not available on this arch\n", papiname);
    return false;
  }

  return true;
}

/*
** is_enabled: determine whether a PAPI counter has already been enabled
**
** Input args: 
**   counter: PAPI counter
**
** Return value: index into papieventlist (success) or negative (not found)
*/
 
int is_enabled (int counter)
{
  int n;

  for (n = 0; n < npapievents; ++n)
    if (papieventlist[n] == counter)
      return n;
  return -1;
}

/*
** enable: enable a PAPI event. ASSUMES that canenable() has already determined
**   that the event can be enabled.
**
** Input args: 
**   counter: PAPI counter
**
** Return value: index into papieventlist
*/
 
int enable (int counter)
{
  int n;

  /* If the event is already enabled, return its index */

  for (n = 0; n < npapievents; ++n) {
    if (papieventlist[n] == counter) {
      return n;
    }
  }

  /* New event */

  papieventlist[npapievents++] = counter;
  return npapievents-1;
}

/*
** getderivedidx: find the table index of a derived counter
**
** Input args: 
**   counter: derived counter
**
** Return value: index into derivedtable (success) or GPTLerror (failure)
*/

int getderivedidx (int dcounter)
{
  int n;

  for (n = 0; n < nderivedentries; ++n) {
    if (derivedtable[n].counter == dcounter)
      return n;
  }
  return GPTLerror ("getderivedidx: failed to find derived counter %d\n", dcounter);
}

/*
** GPTL_PAPIlibraryinit: Call PAPI_library_init if necessary
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTL_PAPIlibraryinit ()
{
  int ret;

  if ((ret = PAPI_is_initialized ()) == PAPI_NOT_INITED) {
    if ((ret = PAPI_library_init (PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
      fprintf(stderr, "GPTL_PAPIlibraryinit: ret=%d PAPI_VER_CURRENT=%d\n", 
	      ret, PAPI_VER_CURRENT);
      return GPTLerror ("GPTL_PAPIlibraryinit: PAPI_library_init failure:%s\n",
			PAPI_strerror (ret));
    }
  }
  return 0;
}

/*
** GPTL_PAPIinitialize(): Initialize the PAPI interface. Called from GPTLinitialize.
**   PAPI_library_init must be called before any other PAPI routines.  
**   PAPI_thread_init is called subsequently if threading is enabled.
**   Finally, allocate space for PAPI counters and start them.
**
** Input args: 
**   maxthreads: number of threads
**
** Return value: 0 (success) or GPTLerror or -1 (failure)
*/
 
int GPTL_PAPIinitialize (const int maxthreads,     /* number of threads */
			 const bool verbose_flag,  /* output verbosity */
			 int *nevents_out,         /* nevents needed by gptl.c */
			 Entry *pr_event_out)      /* events needed by gptl.c */
{
  int ret;       /* return code */
  int n;         /* loop index */
  int t;         /* thread index */
  int *rc;       /* array of return codes from create_and_start_events */
  bool badret;   /* true if any bad return codes were found */

  verbose = verbose_flag;

  /* 
  ** Ensure that PAPI_library_init has already been called.
  */

  if ((ret = GPTL_PAPIlibraryinit ()) < 0)
    return GPTLerror ("GPTL_PAPIinitialize: GPTL_PAPIlibraryinit failure\n");

  /* PAPI_thread_init needs to be called if threading enabled */

#if ( defined THREADED_OMP )
  if (PAPI_thread_init ((unsigned long (*)(void)) (omp_get_thread_num)) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIinitialize: PAPI_thread_init failure\n");
#elif ( defined THREADED_PTHREADS )
  if (PAPI_thread_init ((unsigned long (*)(void)) (pthread_self)) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIinitialize: PAPI_thread_init failure\n");
#endif

  /* allocate and initialize static local space */

  EventSet     = (int *)        GPTLallocate (maxthreads * sizeof (int));
  papicounters = (long_long **) GPTLallocate (maxthreads * sizeof (long_long *));

  for (t = 0; t < maxthreads; t++) {
    EventSet[t] = PAPI_NULL;
    papicounters[t] = (long_long *) GPTLallocate (MAX_AUX * sizeof (long_long));
  }

  /* Event starting apparently must be within a threaded loop. */

  if (npapievents > 0) {
    rc = (int *) GPTLallocate (maxthreads * sizeof (int));

#pragma omp parallel for private (t)

    for (t = 0; t < maxthreads; t++)
      rc[t] = create_and_start_events (t);
  
    badret = false;
    for (t = 0; t < maxthreads; t++)
      if (rc[t] < 0)
	badret = true;
    
    free (rc);

    if (badret)
      return -1;
  }

  *nevents_out = nevents;
  for (n = 0; n < nevents; ++n) {
    pr_event_out[n].counter    = pr_event[n].event.counter;
    pr_event_out[n].counterstr = pr_event[n].event.counterstr;
    pr_event_out[n].prstr      = pr_event[n].event.prstr;
    pr_event_out[n].str        = pr_event[n].event.str;
  }
  return 0;
}

/*
** create_and_start_events: Create and start the PAPI eventset.  File-local.
**   Threaded routine to create the "event set" (PAPI terminology) and start
**   the counters. This is only done once, and is called from GPTL_PAPIinitialize 
** 
** Input args: 
**   t: thread number
**
** Return value: 0 (success) or GPTLerror (failure)
*/

static int create_and_start_events (const int t)  /* thread number */
{
  int ret;
  int n;

  /* Create the event set */

  if ((ret = PAPI_create_eventset (&EventSet[t])) != PAPI_OK)
    return GPTLerror ("create_and_start_events: failure creating eventset: %s\n", 
		      PAPI_strerror (ret));

  /* Add requested events to the event set */

  for (n = 0; n < npapievents; n++) {
    if ((ret = PAPI_add_event (EventSet[t], papieventlist[n])) != PAPI_OK) {
      if (verbose) {
	fprintf (stderr, "%s\n", PAPI_strerror (ret));
	ret = PAPI_event_code_to_name (papieventlist[n], papiname);
	fprintf (stderr, "create_and_start_events: failure adding event:%s\n",
		 papiname);
      }

      if (enable_multiplexing) {
        if (verbose)
	  printf ("Trying multiplexing...\n");
	is_multiplexed = true;
	break;
      } else
	return GPTLerror ("enable_multiplexing is false: giving up\n");
    }
  }

  if (is_multiplexed) {

    /* Cleanup the eventset for multiplexing */

    if ((ret = PAPI_cleanup_eventset (EventSet[t])) != PAPI_OK)
      return GPTLerror ("create_and_start_events: %s\n", PAPI_strerror (ret));
    
    if ((ret = PAPI_destroy_eventset (&EventSet[t])) != PAPI_OK)
      return GPTLerror ("create_and_start_events: %s\n", PAPI_strerror (ret));

    if ((ret = PAPI_create_eventset (&EventSet[t])) != PAPI_OK)
      return GPTLerror ("create_and_start_events: failure creating eventset: %s\n", 
			PAPI_strerror (ret));

    if ((ret = PAPI_multiplex_init ()) != PAPI_OK)
      return GPTLerror ("create_and_start_events: failure from PAPI_multiplex_init%s\n", 
			PAPI_strerror (ret));

    if ((ret = PAPI_set_multiplex (EventSet[t])) != PAPI_OK)
      return GPTLerror ("create_and_start_events: failure from PAPI_set_multiplex: %s\n", 
			PAPI_strerror (ret));

    for (n = 0; n < npapievents; n++) {
      if ((ret = PAPI_add_event (EventSet[t], papieventlist[n])) != PAPI_OK) {
	ret = PAPI_event_code_to_name (papieventlist[n], papiname);
	return GPTLerror ("create_and_start_events: failure adding event:%s\n"
			  "  Error was: %s\n", papiname, PAPI_strerror (ret));
      }
    }
  }

  /* Start the event set.  It will only be read from now on--never stopped */

  if ((ret = PAPI_start (EventSet[t])) != PAPI_OK)
    return GPTLerror ("create_and_start_events: failed to start event set: %s\n", PAPI_strerror (ret));

  return 0;
}

/*
** GPTL_PAPIstart: Start the PAPI counters (actually they are just read).  
**   Called from GPTLstart.
**
** Input args:  
**   t: thread number
**
** Output args: 
**   aux: struct containing the counters
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTL_PAPIstart (const int t,          /* thread number */
		    Papistats *aux)       /* struct containing PAPI stats */
{
  int ret;  /* return code from PAPI lib calls */
  int n;    /* loop index */
  
  /* If no events are to be counted just return */

  if (npapievents == 0)
    return 0;

  /* Read the counters */

  if ((ret = PAPI_read (EventSet[t], papicounters[t])) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIstart: %s\n", PAPI_strerror (ret));

  /* 
  ** Store the counter values.  When GPTL_PAPIstop is called, the counters
  ** will again be read, and differenced with the values saved here.
  */

  for (n = 0; n < npapievents; n++)
    aux->last[n] = papicounters[t][n];
  
  return 0;
}

/*
** GPTL_PAPIstop: Stop the PAPI counters (actually they are just read).  
**   Called from GPTLstop.
**
** Input args:
**   t: thread number
**
** Input/output args: 
**   aux: struct containing the counters
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTL_PAPIstop (const int t,         /* thread number */
		   Papistats *aux)      /* struct containing PAPI stats */
{
  int ret;          /* return code from PAPI lib calls */
  int n;            /* loop index */
  long_long delta;  /* change in counters from previous read */

  /* If no events are to be counted just return */

  if (npapievents == 0)
    return 0;

  /* Read the counters */

  if ((ret = PAPI_read (EventSet[t], papicounters[t])) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIstop: %s\n", PAPI_strerror (ret));
  
  /* 
  ** Accumulate the difference since timer start in aux.
  ** Negative accumulation can happen when multiplexing is enabled, so don't
  ** set count to BADCOUNT in that case.
  */

  for (n = 0; n < npapievents; n++) {
    delta = papicounters[t][n] - aux->last[n];
    if ( ! is_multiplexed && delta < 0)
      aux->accum[n] = BADCOUNT;
    else
      aux->accum[n] += delta;
  }
  return 0;
}

/*
** GPTL_PAPIprstr: Print the descriptive string for all enabled PAPI events.
**   Called from GPTLpr.
**
** Input args: 
**   fp: file descriptor
*/

void GPTL_PAPIprstr (FILE *fp)
{
  int n;
  
  if (narrowprint) {
    for (n = 0; n < nevents; n++) {
      if (strncmp (pr_event[n].event.counterstr, "PAPI_", 5) == 0)
	fprintf (fp, " %8.8s ", &pr_event[n].event.counterstr[5]); /* 5 => lop off "PAPI_" */
      else
	fprintf (fp, " %8.8s ", &pr_event[n].event.counterstr[0]);

      /* Test on < 0 says it's a PAPI counter not a derived counter */

      if (persec && pr_event[n].event.counter < 0)
	fprintf (fp, "  e6/sec");
    }
  } else {
    for (n = 0; n < nevents; n++) {
      fprintf (fp, " %16.16s ", pr_event[n].event.prstr);

      /* Test on < 0 says it's a PAPI counter not a derived counter */

      if (persec && pr_event[n].event.counter < 0)
	fprintf (fp, "  e6/sec");
    }
  }
}

/*
** GPTL_PAPIpr: Print PAPI counter values for all enabled events, including
**   derived events. Called from GPTLpr.
**
** Input args: 
**   fp: file descriptor
**   aux: struct containing the counters
*/

void GPTL_PAPIpr (FILE *fp,                          /* file descriptor to write to */
		  const Papistats *aux,              /* stats to write */
		  const int t,                       /* thread number */
		  const int count,                   /* number of invocations */
		  const double wcsec)                /* wallclock time (sec) */
{
  const char *shortintfmt   = "%8ld ";
  const char *longintfmt    = "%16ld ";
  const char *shortfloatfmt = "%8.2e ";
  const char *longfloatfmt  = "%16.10e ";
  const char *intfmt;       /* integer format */
  const char *floatfmt;     /* floating point format */

  int n;              /* loop index */
  int numidx;         /* index pointer to appropriated (derived) numerator */
  int denomidx;       /* index pointer to appropriated (derived) denominator */
  double val;         /* value to be printed */

  intfmt   = narrowprint ? shortintfmt   : longintfmt;
  floatfmt = narrowprint ? shortfloatfmt : longfloatfmt;

  for (n = 0; n < nevents; n++) {
    numidx = pr_event[n].numidx;
    if (pr_event[n].event.counter > 0) {   /* derived event */
      denomidx = pr_event[n].denomidx;

      /* Protect against divide by zero */

      if (aux->accum[denomidx] > 0)
	val = (double) aux->accum[numidx] / (double) aux->accum[denomidx];
      else
	val = 0.;
      fprintf (fp, floatfmt, val);

    } else {                               /* Raw PAPI event */

      if (aux->accum[numidx] < PRTHRESH)
	fprintf (fp, intfmt, (long) aux->accum[numidx]);
      else
	fprintf (fp, floatfmt, (double) aux->accum[numidx]);

      if (persec) {
	if (wcsec > 0.)
	  fprintf (fp, "%8.2f ", aux->accum[numidx] * 1.e-6 / wcsec);
	else
	  fprintf (fp, "%8.2f ", 0.);
      }
    }
  }
}

/*
** GPTL_PAPIprintenabled: Print list of enabled timers
**
** Input args:
**   fp: file descriptor
*/

void GPTL_PAPIprintenabled (FILE *fp)
{
  int n;

  fprintf (fp, "PAPI events enabled (including those required for derived events):\n");
  for (n = 0; n < nevents; n++)
    fprintf (fp, "  %s\n", pr_event[n].event.str);
  fprintf (fp, "\n");
}  

/*
** GPTL_PAPIadd: Accumulate PAPI counters. Called from add.
**
** Input/Output args: 
**   auxout: auxout = auxout + auxin
**
** Input args:
**   auxin: counters to be summed into auxout
*/

void GPTL_PAPIadd (Papistats *auxout,      /* output struct */
		   const Papistats *auxin) /* input struct */
{
  int n;
  
  for (n = 0; n < npapievents; n++)
    if (auxin->accum[n] == BADCOUNT || auxout->accum[n] == BADCOUNT)
      auxout->accum[n] = BADCOUNT;
    else
      auxout->accum[n] += auxin->accum[n];
}

/*
** GPTL_PAPIfinalize: finalization routine must be called from single-threaded
**   region. Free all malloc'd space
*/

void GPTL_PAPIfinalize (int maxthreads)
{
  int t;

  for (t = 0; t < maxthreads; t++) {
    free (papicounters[t]);
  }

  free (EventSet);
  free (papicounters);

  /* Reset initial values */

  nevents = 0;
  npapievents = 0;
}

/*
** GPTL_PAPIquery: return current PAPI counter info. Return into a long for best
**   compatibility possibilities with Fortran.
**
** Input args:
**   aux:       struct containing the counters
**   ncounters: max number of counters to return
**
** Output args:
**   papicounters_out: current value of PAPI counters
*/

void GPTL_PAPIquery (const Papistats *aux,
		     long long *papicounters_out,
		     int ncounters)
{
  int n;

  if (ncounters > 0) {
    for (n = 0; n < ncounters && n < npapievents; n++) {
      papicounters_out[n] = (long long) aux->accum[n];
    }
  }
}

/*
** GPTL_PAPIis_multiplexed: return status of whether events are being multiplexed
*/

bool GPTL_PAPIis_multiplexed ()
{
  return is_multiplexed;
}

/*
** The following functions are publicly available
*/

void read_counters100 ()
{
  int i;
  int ret;
  long_long counters[MAX_AUX];

  for (i = 0; i < 10; ++i) {
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
    ret = PAPI_read (EventSet[0], counters);
  }
  return;
}

/*
** GPTL_PAPIname2id: convert a PAPI event name in string form to an int
**
** Input args:
**   name: PAPI event name
**
** Return value: event index (success) or 77 (error)
*/

int GPTL_PAPIname2id (const char *name)
{
  int counter, n, ret;

  /* Check for a GPTL derived event */
  for (n = 0; n < nderivedentries; n++) {
  if (strcmp (name, derivedtable[n].counterstr) == 0)
    return derivedtable[n].counter;
  }

  /* If not, check PAPI presets */
  for (n = 0; n < npapientries; n++) {
  if (strcmp (name, papitable[n].counterstr) == 0)
    return papitable[n].counter;
  }

  /* If not, ... */
  /* Initialize PAPI if it hasn't already been done */
  if (GPTL_PAPIlibraryinit () < 0)
    return GPTLerror ("GPTL_PAPIsetoption: PAPI library init error\n");

  /* Check native PAPI events */
  if ((ret = PAPI_event_name_to_code (name, &counter)) != PAPI_OK)
    return GPTLerror ("GPTL_PAPIname2id: PAPI_strerror: %s\n", PAPI_strerror (ret));

  return counter; 

}
#else

/*
** "Should not be called" entry points for publicly available GPTL_PAPI routines
*/

#include <stdio.h>

/*
** GPTL_PAPIlibraryinit: Call PAPI_library_init if necessary
**
** Return value: 0 (success) or GPTLerror (failure)
*/

void GPTL_PAPIlibraryinit ()
{
  printf ("PAPI not enabled: GPTL_PAPIlibraryinit should not be called\n");
}

int GPTL_PAPIname2id (const char *name, int nc)
{
  printf ("PAPI not enabled: GPTL_PAPIname2id should not be called\n");
  return 77;
}

#endif  /* HAVE_PAPI */

