#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

TEST_CASE("pass", "[fake_infra_test]")
{
#ifdef SCREAM_FORCE_RUN_FAIL
  REQUIRE(false);
#else
  REQUIRE(true);
#endif
}

} // empty namespace
