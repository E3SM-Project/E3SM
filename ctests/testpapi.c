#include "../gptl.h"
#include <stdio.h>

int main (int argc, char **argv)
{
  int ret;
  int i, code;
  long long pc[1]; /* papi counters */
  double sum;

  printf ("testpapi: Testing PAPI interface...\n");

  printf ("%s: testing getting event code for PAPI_TOT_CYC...\n", argv[0]);
  if ((ret = GPTLevent_name_to_code ("PAPI_TOT_CYC", &code)) != 0) {
    printf ("Failure\n");
    return 2;
  }
  printf ("Success\n");

  printf ("%s: testing GPTLsetoption(PAPI_TOT_CYC,1)...\n", argv[0]);
  if (GPTLsetoption (code, 1) != 0) {
    printf ("Failure\n");
    return 3;
  }
  printf ("Success\n");

  printf ("%s: testing GPTLinitialize\n", argv[0]);
  if ((ret = GPTLinitialize ()) != 0) {
    printf ("Failure\n");
    return 3;
  }
  printf ("Success\n");

  printf ("%s: testing GPTLstart\n", argv[0]);
  if ((ret = GPTLstart ("sum")) != 0) {
    printf ("Unexpected failure from GPTLstart(\"sum\")\n");
    return 3;
  }
  printf ("Success\n");

  sum = 0.;
  for (i = 0; i < 1000000; ++i) 
    sum += (double) i;
  printf ("%s: testing GPTLstop\n", argv[0]);
  if ((ret = GPTLstop ("sum")) != 0) {
    printf ("Unexpected failure from GPTLstop(\"sum\")\n");
    return 3;
  }
  printf ("Success\n");

  printf ("%s: testing GPTLquerycounters...\n", argv[0]);
  if (GPTLquerycounters ("sum", 0, pc) != 0) {
    printf ("Failure\n");
    return 4;
  }
  printf ("sum=%g\n",sum);
  printf ("%s: testing reasonableness of PAPI counters...\n", argv[0]);
  if (pc[0] < 1 || pc[0] > 1.e8) {
    printf ("Suspicious PAPI_TOT_CYC value=%lld for 1e6 additions\n", pc[0]);
    return 5;
  } else {
    printf ("Success\n");
  }
  printf ("%s: all tests successful\n", argv[0]);
  return 0;
}
