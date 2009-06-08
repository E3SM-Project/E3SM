/*
** print_memusage:
**
**   Prints info about memory usage of this process by calling get_memusage.
**
**   Return value: 0  = success
**                 -1 = failure
*/

#include <stdio.h>
#include "gptl.h"

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
