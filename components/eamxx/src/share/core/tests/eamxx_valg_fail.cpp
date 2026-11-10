#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

TEST_CASE("force_valgrind_err")
{
  bool* uninit = new bool[1];
  int i = 0;
  if (uninit[0]) {
    ++i;
  }
  else {
    i += 4;
  }
  if (i<4) {
    printf("less than four\n");
  }
  delete uninit;
}

} // empty namespace
