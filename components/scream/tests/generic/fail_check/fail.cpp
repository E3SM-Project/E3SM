#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

TEST_CASE("force_fail")
{
  REQUIRE(false); // force this test to fail
}

} // empty namespace
