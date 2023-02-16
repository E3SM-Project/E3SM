#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

TEST_CASE("force_valgrind_err")
{
  bool uninit;
  int i = 0;
  if (uninit) {
    ++i;
  }
  else {
    i += 4;
  }
  REQUIRE(i < 10);
}

} // empty namespace
