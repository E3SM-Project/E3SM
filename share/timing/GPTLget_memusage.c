/*
** $Id: get_memusage.c,v 1.10 2010-11-09 19:08:53 rosinski Exp $
**
** Author: Jim Rosinski
**   Credit to Chuck Bardeen for MACOS section (__APPLE__ ifdef)
**
** get_memusage:
**
**   Designed to be called from Fortran, returns information about memory
**   usage in each of 5 input long long* args.  On Linux read from the /proc
**   filesystem because getrusage() returns placebos (zeros).  Return -1 for
**   values which are unavailable or ambiguous on a particular architecture.
**   Reported numbers are in kilobytes.
**
**   Return value: 0  = success
**                 -1 = failure
*/

#include <sys/resource.h>
#include "gptl.h"    /* additional cpp defs and function prototypes */
#include <stdio.h>

/* _AIX is automatically defined when using the AIX C compilers */
#ifdef _AIX
#include <sys/times.h>
#endif

#ifdef IRIX64
#include <sys/time.h>
#endif

#ifdef HAVE_SLASHPROC

#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#elif (defined __APPLE__)

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#endif

#ifdef BGP

#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <malloc.h>
#define   Personality                    _BGP_Personality_t

#endif

#ifdef __bgq__

#include <malloc.h>
#include <spi/include/kernel/memory.h>

#endif

#define PRINT_MEMUSAGE 0

int GPTLget_memusage (long long *size, long long *rss, long long *share, long long *text, long long *datastack)
{
#if defined (BGP)
  long long alloc, total;
  int node_config;
  struct mallinfo m;
  Personality pers;

 /* memory available */
  Kernel_GetPersonality(&pers, sizeof(pers));
  total = BGP_Personality_DDRSizeMB(&pers);

  node_config  = BGP_Personality_processConfig(&pers);
  if (node_config == _BGP_PERS_PROCESSCONFIG_VNM) total /= 4;
  else if (node_config == _BGP_PERS_PROCESSCONFIG_2x2) total /= 2;
  total *= 1024; // in KB

  /* total memory used  - heap only (not static memory)*/
  *size = total;

  m = mallinfo();
  alloc = m.hblkhd + m.uordblks;

  *rss = alloc;
  *share     = -1;
  *text     = -1;
  *datastack = -1;
  return 0;

#elif (defined __bgq__)
  uint64_t heap, shared, stack;

  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);

  *size      = heap/1024;
  *rss       = heap/1024;
  *share     = shared/1024;
  *text      = -1;
  *datastack = stack/1024;
  return 0;

#elif (defined HAVE_SLASHPROC)
  FILE *fd;                       /* file descriptor for fopen */
  int pid;                        /* process id */
  char file[32];                  /* full path to file in /proc */
  char line[128];                 /* line buffer for reading /proc/status */

  /*
  ** Read from /proc/<pid>/status which reports labeled fields in KB.
  ** Use VmHWM (peak RSS) for size instead of VmSize (virtual address space)
  ** because VmSize is inflated by GPU BAR mappings and other non-physical
  ** reservations, while VmHWM is the true memory high-water mark.
  */

  pid = (int) getpid ();
  if (pid <= 0) {
    fprintf (stderr, "get_memusage: pid %d is non-positive\n", pid);
    return -1;
  }

  *size      = -1;
  *rss       = -1;
  *share     = -1;
  *text      = -1;
  *datastack = -1;

  sprintf (file, "/proc/%d/status", pid);
  if ((fd = fopen (file, "r")) == NULL) {
    fprintf (stderr, "get_memusage: bad attempt to open %s\n", file);
    return -1;
  }

  while (fgets (line, sizeof(line), fd)) {
    if      (sscanf (line, "VmHWM: %lld", size)      == 1) ;
    else if (sscanf (line, "VmRSS: %lld", rss)       == 1) ;
    else if (sscanf (line, "VmData: %lld", datastack) == 1) ;
    else if (sscanf (line, "VmExe: %lld", text)       == 1) ;
    else if (sscanf (line, "RssFile: %lld", share)    == 1) ;
  }
  fclose (fd);

#if PRINT_MEMUSAGE
  fprintf (stderr, "get_memusage: size=%lld KB, rss=%lld KB, share=%lld KB, text=%lld KB, datastack=%lld KB\n",
           *size, *rss, *share, *text, *datastack);
#endif

  return 0;

#elif (defined __APPLE__)

  FILE *fd;
  char cmd[60];
  int pid = (int) getpid ();

  // returned values are in KBs
  sprintf (cmd, "ps -o vsz -o rss -o tsiz -p %d | grep -v RSS", pid);
  fd = popen (cmd, "r");

  if (fd) {
    fscanf (fd, "%lld %lld %lld", size, rss, text);
    *share     = -1;
    *datastack = -1;
    (void) pclose (fd);
  }

  return 0;

#else

  struct rusage usage;         /* structure filled in by getrusage */

  if (getrusage (RUSAGE_SELF, &usage) < 0)
    return -1;

  *size      = -1;
  *rss       = (long long) usage.ru_maxrss; // in KBs
  *share     = -1;
  *text      = -1;
  *datastack = -1;
#ifdef IRIX64
  *datastack = (long long) (usage.ru_idrss + usage.ru_isrss);
#endif
  return 0;

#endif
}
