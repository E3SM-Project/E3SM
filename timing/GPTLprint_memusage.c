/*
** $Id: print_memusage.c,v 1.5 2009/03/24 20:51:25 rosinski Exp $
**
** Author: Jim Rosinski
**
** print_memusage:
**
**   Prints info about memory usage of this process by calling get_memusage.
**
**   Return value: 0  = success
**                 -1 = failure
*/

#include "gptl.h"
#include <stdio.h>

int GPTLprint_memusage (const char *str)
{
  int size;
  int rss;
  int share;
  int text;
  int datastack;

  if (GPTLget_memusage (&size, &rss, &share, &text, &datastack) < 0)
    return -1;

  printf ("%s size=%d rss=%d share=%d text=%d datastack=%d\n", 
	  str, size, rss, share, text, datastack);
  return 0;
}
