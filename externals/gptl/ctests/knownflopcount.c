#include <unistd.h>  /* getopt */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>

int main (int argc, char **argv)
{
  const int niter = 1000;            /* iteration count */
  const int arrlen = 1000000;        /* array size */
  const int maxevents = 10;          /* Max PAPI events */
  long_long *papicounters;           /* output from PAPI_read */
  long_long *prvcounters;            /* output from PAPI_read */
  long_long diff;                    /* PAPI count for region */
  char eventname[maxevents][PAPI_MAX_STR_LEN];
  
  int EventSet = PAPI_NULL;  /* Event set needed by PAPI lib */

  int ret;            /* return code */
  int code;           /* PAPI event code */
  int i, n;           /* loop indices */
  int c;              /* for parsing argv */
  int nevents = 0;    /* number of PAPI events (init to 0) */
    
  double arr[arrlen]; /* array to do math on */

  void init (double *);  /* initialize arr */

  papicounters = (long_long *) malloc (maxevents * sizeof (long_long));
  prvcounters  = (long_long *) malloc (maxevents * sizeof (long_long));

  /* Initialize the PAPI library */

  if ((ret = PAPI_library_init (PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
    printf ("%s\n", PAPI_strerror (ret));
    return -1;
  }

  /* Create the eventset */

  if ((ret = PAPI_create_eventset (&EventSet)) != PAPI_OK) {
    printf ("Failure creating eventset: %s\n", PAPI_strerror (ret));
    return -1;
  }

  if (argc < 2) {
    printf ("Usage: %s -e papi_counter_name ...\n", argv[0]);
    return -1;
  }

  /* Parse arg list */

  while ((c = getopt (argc, argv, "e:")) != -1) {
    switch (c) {
    case 'e':

      /* Convert name to code */

      if ((ret = PAPI_event_name_to_code (optarg, &code)) != PAPI_OK) {
	printf ("No code found for event %s\n", optarg);
	printf ("PAPI_strerror says: %s\n", PAPI_strerror (ret));
	return -1;
      }

      /* Add the event */

      if ((ret = PAPI_add_event (EventSet, code)) != PAPI_OK) {
	printf ("%s\n", PAPI_strerror (ret));
	printf ("Failure adding event %s\n", optarg);
	return -1;
      }

      if (nevents >= maxevents) {
	printf ("%d events is too many\n", nevents);
	return -1;
      }
      strcpy (eventname[nevents++], optarg);
      break;
    default:
      printf ("unknown option %c\n", c);
      return -1;
    }
  }

  /* Start the eventset */

  if ((ret = PAPI_start (EventSet)) != PAPI_OK)
    printf ("%s\n", PAPI_strerror (ret));

  init (arr);

  /* Read counters before computation */

  if ((ret = PAPI_read (EventSet, prvcounters)) != PAPI_OK) {
    printf ("PAPI_read error\n");
    return -1;
  }

  /* Do computation */

  for (n = 0; n < niter; ++n)
    for (i = 0; i < arrlen; ++i)
      arr[i] += 0.1e0*arr[i];

  /* Read counters after computation */

  if ((ret = PAPI_read (EventSet, papicounters)) != PAPI_OK) {
    printf ("PAPI_read error\n");
    return -1;
  }

  /* Print counter information */

  printf ("FP_OPS and FP_INS should be 2.e9\n\n");

  for (n = 0; n < nevents; ++n) {
    diff = papicounters[n] - prvcounters[n];
    printf ("%s count = %20.14e\n", eventname[n], (double) diff);
  }
  return 0;
}

void init (double *arr)
{
  int i;

  for (i = 0; i < 1000000; ++i)
    arr[i] = 1./i;
}

