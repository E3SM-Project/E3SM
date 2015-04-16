/*
** print_rusage.c
**
** Author: Jim Rosinski
**
** print_rusage:
**
**   Prints info from getrusage()
**
**   Return value: 0  = success
**                 -1 = failure
*/

#include "private.h"
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>

int GPTLprint_rusage (const char *str)
{
  struct rusage usage;
  static const char *thisfunc = "GPTLprint_rusage";
  static const float onek = 1042.;
  
  if (getrusage (RUSAGE_SELF, &usage) < 0)
    return GPTLerror ("%s: failure from getrusage()\n", thisfunc);

  /* ru_maxrss is in KB */
  printf ("%s ru_maxrss=%.1f MB ru_minflt=%.1f K ru_majflt=%.1f K ru_nvcsw=%.1f K\n", str,
	  usage.ru_maxrss/onek, 
	  usage.ru_minflt/onek, 
	  usage.ru_majflt/onek, 
	  usage.ru_nvcsw/onek);
  return 0;
}
