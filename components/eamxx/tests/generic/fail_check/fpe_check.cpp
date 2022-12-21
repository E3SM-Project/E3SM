#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

TEST_CASE("force_fpe")
{
  float foo = 42.0;
  float bar = foo / 0.0;
  std::cout << bar << std::endl;

  REQUIRE(true);
}

} // empty namespace
