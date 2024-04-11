#include <iostream>

int main(int, char**)
{
#ifdef SCREAM_FORCE_RUN_FAIL
  return 1;
#else
  return 0;
#endif
}
