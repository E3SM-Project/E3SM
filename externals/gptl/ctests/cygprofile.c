#include <stdio.h>
#include <stdlib.h>
#include "../gptl.h"

extern void callsubs (int);

int main (int argc, char **argv)
{
  int niter = 1000;
  int ret;
  int n, nregions;
  char name[64];
  double wallclock;
  
  if (argc == 2) {
    niter = atoi (argv[1]);
  } else if (argc > 2) {
    printf ("Usage: %s loop_length\n", argv[0]);
  }

  GPTLsetoption (GPTLabort_on_error, 0);

  if ((ret = GPTLinitialize ()) != 0) {
    printf ("%s: GPTLinitialize failure\n", argv[0]);
    return -1;
  }
  callsubs (niter);
  GPTLpr (0);

  printf ("%s: Testing GPTLget_nregions...\n", argv[0]);
  if (GPTLget_nregions (0, &nregions) < 0) {
    printf ("%s: GPTLget_nregions failure\n", argv[0]);
    return -1;
  }
  printf ("Success\n");

  printf ("%s: Testing GPTLget_regionname, GPTLget_wallclock "
	  "for %d regions...\n", argv[0], nregions);
  for (n = 0; n < nregions; ++n) {
    if ((ret = GPTLget_regionname (0, n, name, sizeof (name))) < 0) {
      printf ("%s: GPTLget_regionname failure\n", argv[0]);
      return -1;
    }

    if ((ret = GPTLget_wallclock (name, 0, &wallclock)) < 0) {
      printf ("%s: GPTLget_wallclock failure for name=%s\n", argv[0], name);
      return -1;
    }
  }
  printf ("Success\n");

  (void) GPTLfinalize ();
  return 0;
}
