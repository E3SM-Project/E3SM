/*
** $Id: print_memusage.c,v 1.13 2010-11-09 19:08:54 rosinski Exp $
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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef __bgq__
#include <spi/include/kernel/memory.h>
#endif

static int convert_to_mb = 1;   /* true */

int GPTLprint_memusage (const char *str)
{
#ifdef __bgq__
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;

  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

  printf("%s Memory(MB): heap-alloc: %.2f, heap-avail: %.2f,"
         "stack-alloc: %.2f, stack-avail: %.2f,"
         "shared: %.2f, persist: %.2f, guard: %.2f, mmap: %.2f\n", str,
         (double)heap/(1024*1024),   (double)heapavail/(1024*1024),
         (double)stack/(1024*1024),  (double)stackavail/(1024*1024),
         (double)shared/(1024*1024), (double)persist/(1024*1024),
         (double)guard/(1024*1024),  (double)mmap/(1024*1024));
  return 0;

#else
  long long size;                           /* peak RSS in KB (VmHWM on Linux) */
  long long rss;                            /* resident set size in KB */
  long long share;                          /* shared pages in KB */
  long long text;                           /* text segment in KB */
  long long datastack;                      /* data segment in KB */

  setbuf(stdout, NULL); // don't buffer stdout, flush
  if (GPTLget_memusage (&size, &rss, &share, &text, &datastack) < 0) {
    printf ("GPTLprint_memusage: GPTLget_memusage failed.\n");
    return -1;
  }

  /* GPTLget_memusage returns values in KB on all platforms */
  if (convert_to_mb)
    printf ("%s sysmem size=%.1f MB rss=%.1f MB share=%.1f MB text=%.1f MB datastack=%.1f MB\n",
            str, size/1024., rss/1024., share/1024., text/1024., datastack/1024.);
  else
    printf ("%s sysmem size=%lld rss=%lld share=%lld text=%lld datastack=%lld\n",
            str, size, rss, share, text, datastack);

  return 0;
#endif
}

