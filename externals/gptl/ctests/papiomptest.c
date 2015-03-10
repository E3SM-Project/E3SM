#include <stdio.h>
#include <stdlib.h>  /* atoi,exit */
#include <unistd.h>  /* getopt */
#include <string.h>  /* memset */
#ifdef THREADED_OMP
#include <omp.h>
#endif

#include "../gptl.h"

double add (int, double);
double multiply (int, int, double);
double multadd (int, double);
double divide (int, double);
double compare (int, int);

int main (int argc, char **argv)
{
  int nompiter = 128;
  int looplen = 1000000;
  int iter;
  int papiopt;
  int c;
  int ret;
  double value;
  extern char *optarg;

  printf ("Purpose: test known-length loops with various floating point ops\n");
  printf ("Include PAPI and OpenMP, respectively, if enabled\n");
  printf ("Usage: %s [-l looplen] [-n nompiter] [-p papi_option_name]\n", argv[0]);

  while ((c = getopt (argc, argv, "l:n:p:")) != -1) {
    switch (c) {
    case 'l':
      looplen = atoi (optarg);
      printf ("Set looplen=%d\n", looplen);
      break;
    case 'n':
      nompiter = atoi (optarg);
      printf ("Set nompiter=%d\n", nompiter);
      break;
    case 'p':
      if ((ret = GPTLevent_name_to_code (optarg, &papiopt)) != 0) {
	printf ("Failure from GPTLevent_name_to_code(%s)\n", optarg);
	exit (1);
      }
      
      if (GPTLsetoption (papiopt, 1) < 0) {
	printf ("Failure from GPTLsetoption (%s,1)\n", optarg);
	exit (1);
      }
      break;
    default:
      printf ("unknown option %c\n", c);
      exit (2);
    }
  }
  
  printf ("Outer loop length (OMP)=%d\n", nompiter);
  printf ("Inner loop length=%d\n", looplen);

  ret = GPTLsetoption (GPTLverbose, 1);
  ret = GPTLsetoption (GPTLabort_on_error, 1);
  ret = GPTLsetoption (GPTLoverhead, 1);
  ret = GPTLsetoption (GPTLnarrowprint, 1);

  if ((GPTLinitialize ()) != 0) {
    printf ("papiomptest: GPTLinitialize failure\n");
    return -1;
  }

  ret = GPTLstart ("total");
#ifdef THREADED_OMP
#pragma omp parallel for private (iter, value)
#endif
  for (iter = 1; iter <= nompiter; iter++) {
    value = add (looplen, 0.);
    value = multiply (looplen, iter, 0.);
    value = multadd (looplen, 0.);
    value = divide (looplen, 1.);
    value = compare (looplen, iter);
  }
  ret = GPTLstop ("total");
  ret = GPTLpr (0);
  if (GPTLfinalize () < 0)
    exit (6);

  return 0;
}

double add (int looplen, double zero)
{
  int i;
  char string[128];
  double val = zero;

  if (looplen < 1000)
    sprintf (string, "%dadditions", looplen);
  else
    sprintf (string, "%-.3gadditions", (double) looplen);

  if (GPTLstart (string) < 0)
    exit (1);

  for (i = 1; i <= looplen; ++i)
      val += zero;

  if (GPTLstop (string) < 0)
    exit (1);

  return val;
}

double multiply (int looplen, int iter, double zero)
{
  int i;
  char string[128];
  double val = iter;

  if (looplen < 1000)
    sprintf (string, "%dmultiplies", looplen);
  else
    sprintf (string, "%-.3gmultiplies", (double) looplen);

  if (GPTLstart (string) < 0)
    exit (1);

  for (i = 1; i <= looplen; ++i) {
    val *= zero;
  }

  if (GPTLstop (string) < 0)
    exit (1);

  return val;
}

double multadd (int looplen, double zero)
{
  int i;
  char string[128];
  double val = zero;

  if (looplen < 1000)
    sprintf (string, "%dmultadds", looplen);
  else
    sprintf (string, "%-.3gmultadds", (double) looplen);

  if (GPTLstart (string) < 0)
    exit (1);

  for (i = 1; i <= looplen; ++i)
    val += zero * i;

  if (GPTLstop (string) < 0)
    exit (1);

  return val;
}

double divide (int looplen, double one)
{
  int i;
  char string[128];
  double val = one;

  if (looplen < 1000)
    sprintf (string, "%ddivides", looplen);
  else
    sprintf (string, "%-.3gdivides", (double) looplen);

  if (GPTLstart (string) < 0)
    exit (1);

  for (i = 1; i <= looplen; ++i)
    val /= one;

  if (GPTLstop (string) < 0)
    exit (1);

  return val;
}

double compare (int looplen, int iter)
{
  int i;
  char string[128];
  double val = iter;

  if (looplen < 1000)
    sprintf (string, "%dcompares", looplen);
  else
    sprintf (string, "%-.3gcompares", (double) looplen);

  if (GPTLstart (string) < 0)
    exit (1);

  for (i = 0; i < looplen; ++i)
    if (val < i)
      val = i;

  if (GPTLstop (string) < 0)
    exit (1);

  return val;
}
