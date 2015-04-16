#include "../gptl.h"
#ifdef HAVE_PAPI
#include <papi.h>
#endif
#include "myclasses.h"

int main ()
{
  X *x;
  Y *y;
  int ret;

#ifdef HAVE_PAPI
  GPTLsetoption(GPTLmultiplex,1);
  GPTLsetoption(PAPI_L2_DCH, 1);
  GPTLsetoption(PAPI_L1_TCM, 1);
  GPTLsetoption(PAPI_L3_TCM, 1);
#endif

  ret = GPTLinitialize ();
  ret = GPTLstart ("total");
  x = new (X);
  x->func (1.2);
  x->func (1);
  delete (x);

  y = new (Y);
  y->func (1.2);
  y->func (1);
  delete (y);

  ret = GPTLstop ("total");
  ret = GPTLpr (0);
}

