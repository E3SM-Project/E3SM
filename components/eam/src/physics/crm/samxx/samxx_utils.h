
#pragma once

#include <type_traits>

#include "samxx_const.h"
#include "vars.h"

extern "C" {
 void scream_session_init();
 void scream_session_finalize();
}

YAKL_INLINE int is_same_str(const char *str_a, const char *str_b)
{
  int match = 0;
  unsigned i = 0;
  unsigned done = 0;
  while ((match == 0) && !done)
  {
    if ((str_a[i] == 0) || (str_b[i] == 0)) {
      done = 1;
    } else if (str_a[i] != str_b[i]) {
      match = i+1;
      if (((int)str_a[i] - (int)str_b[i]) < 0) match = 0 - (i + 1);
    }
    i++;
  }
  return match;
}

