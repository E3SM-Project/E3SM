#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

#ifdef SCREAM_FORCE_BUILD_FAIL
#error "Forcing failure to test test-all-scream"
#endif

#ifdef SCREAM_FORCE_RUN_FAIL
TEST_CASE("force_fail", "[fake_infra_test]")
{
  REQUIRE(false); // force this test to fail
}
#endif

TEST_CASE("pass", "[fake_infra_test]")
{
  REQUIRE(true);
}

} // empty namespace
