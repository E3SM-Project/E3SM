#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

#ifdef SCREAM_FORCE_BUILD_FAIL
#error "Forcing failure to test test-all-scream"
#endif

TEST_CASE("pass", "[fake_infra_test]")
{
  REQUIRE(true);
}

} // empty namespace
