/*
** get_memusage: 
**
**   Designed to be called from Fortran, returns information about memory
**   usage in each of 5 input int* args.  On Linux read from the /proc
**   filesystem because getrusage() returns placebos (zeros).  Return -1 for
**   values which are unavailable or ambiguous on a particular architecture.
**
**   Return value: 0  = success
**                 -1 = failure
*/

#include <sys/resource.h>

#ifdef AIX
#include <sys/times.h>
#endif

#ifdef IRIX64
#include <sys/time.h>
#endif

#ifdef LINUX
#include <sys/time.h>
#include <sys/types.h>
#include <stdio.h>
#include <unistd.h>
#endif

int GPTLget_memusage (int *size, int *rss, int *share, int *text, int *datastack)
{
#ifdef LINUX
  FILE *fd;                       /* file descriptor for fopen */
  int pid;                        /* process id */
  static char *head = "/proc/";   /* part of path */
  static char *tail = "/statm";   /* part of path */
  char file[19];                  /* full path to file in /proc */
  int dum;                        /* unused */

  /*
  ** The file we want to open is /proc/<pid>/statm
  */

  pid = (int) getpid ();
  if (pid > 999999) {
    fprintf (stderr, "get_memusage: pid %d is too large\n", pid);
    return -1;
  }

  sprintf (file, "%s%d%s", head, pid, tail);
  if ((fd = fopen (file, "r")) < 0) {
    fprintf (stderr, "get_memusage: bad attempt to open %s\n", file);
    return -1;
  }

  /*
  ** Read the desired data from the /proc filesystem directly into the output
  ** arguments, close the file and return.
  */

  (void) fscanf (fd, "%d %d %d %d %d %d %d", size, rss, share, text, datastack, &dum, &dum);
  (void) fclose (fd);
  return 0;

#else

  struct rusage usage;         /* structure filled in by getrusage */

  if (getrusage (RUSAGE_SELF, &usage) < 0)
    return -1;
  
  *size      = -1;
  *rss       = usage.ru_maxrss;
  *share     = -1;
  *text      = -1;
  *datastack = -1;
#ifdef IRIX64
  *datastack = usage.ru_idrss + usage.ru_isrss;
#endif
  return 0;

#endif
}
