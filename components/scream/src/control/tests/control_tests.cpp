#include <catch2/catch.hpp>
#include <control/atmosphere_driver.hpp>

namespace {

TEST_CASE("Atmosphere Driver", "stub") {
  using namespace scream::control;

  int val = driver_stub();
  REQUIRE(val == 42);
}

} // empty namespace
