#include <catch2/catch.hpp>

TEST_CASE("force_fail")
{
  REQUIRE(false); // force this test to fail
}
