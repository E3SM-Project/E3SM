#include "../gptl.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include <stdio.h>

int handle; /* for _handle routines--used by both granularity and overhead */

int main (int argc, char **argv)
{
  typedef struct {
    char *name;
    int utr;
  } Vals;
  /*
  ** Don't test MPI_Wtime: too complex to get to work everywhere.
  */
  Vals vals[] = {{"gettimeofday",   GPTLgettimeofday},
		 {"nanotime",       GPTLnanotime},
		 /*		 {"mpiwtime",       GPTLmpiwtime}, */
		 {"clockgettime",   GPTLclockgettime},
		 {"papitime",       GPTLpapitime},
		 {"read_real_time", GPTLread_real_time}};
  static const int nvals = sizeof (vals) / sizeof (Vals);
    
  int ret;
  int n;
  char *mingranname = 0;
  char *minohname = 0;
  double mingran = 99999.;
  double minoh = 99999.;
  double val;

  extern double granularity (void);
  extern double overhead (void);

  for (n = 0; n < nvals; n++) {
    printf ("Checking %s...\n", vals[n].name);
    if ((ret = GPTLsetutr (vals[n].utr)) == 0) {
      ret = GPTLsetoption (GPTLoverhead, 0);

      if ((ret = GPTLinitialize ()) != 0) {
	printf ("GPTLinitialize failure\n");
	return -1;
      }

      /* Need to call MPI_Init when the UTR is MPI_Wtime */

      /*
      ** Don't test MPI_Wtime: too complex to get to work everywhere.
      */

      /*
      if (vals[n].utr == GPTLmpiwtime) {
	if ((ret = MPI_Init (&argc, &argv)) != 0) {
	  printf ("Failure from MPI_Init: skipping MPI_Wtime...\n");
	  continue;
	}
      }
      */

      /*
      ** Warm up the handle routines, especially the first "start" call
      ** which adds the entry point.
      */

      handle = 0;
      ret = GPTLstart_handle ("zzz", &handle);
      ret = GPTLstop_handle ("zzz", &handle);
      ret = GPTLreset ();

      val = granularity ();
      if (val < mingran) {
	mingran = val;
	mingranname = vals[n].name;
      }
	
      ret = GPTLreset ();
      val = overhead ();
      if (val < minoh) {
	minoh = val;
	minohname = vals[n].name;
      }
    } else {
      printf ("Not available\n");
    }
    ret = GPTLfinalize ();
    printf ("\n");
  }
  printf ("func with finest granularity = %s (%g)\n", mingranname, mingran);
  printf ("func with min overhead       = %s (%g)\n", minohname, minoh);
  return 0;
}

double granularity ()
{
  int ret;
  int count;
  int onflg;
  double wallclock;
  double usr;
  double sys;
  long long papi;

  /* handle was initialized in main*/

  do {
    ret = GPTLstart_handle ("zzz", &handle);
    ret = GPTLstop_handle ("zzz", &handle);
    ret = GPTLquery ("zzz", 0, &count, &onflg, &wallclock, &usr, &sys, &papi, 0);
  } while (wallclock == 0.);
  
  printf ("granularity = %g seconds found after %d iterations\n", wallclock, count);
  return wallclock;
}

double overhead ()
{
  int n;
  double oh;

  int ret;
  int count;
  int onflg;
  double wallclock;
  double usr;
  double sys;
  long long papi;

  /* handle was initialized in main*/

  for (n = 0; n < 10000; n++) {
    ret = GPTLstart_handle ("zzz", &handle);
    ret = GPTLstop_handle ("zzz", &handle);
  }
  ret = GPTLquery ("zzz", 0, &count, &onflg, &wallclock, &usr, &sys, &papi, 0);
  oh = 0.0001 * wallclock;
  printf ("overhead = %g per call based on 10,000 iterations\n", oh);
  return oh;
}
