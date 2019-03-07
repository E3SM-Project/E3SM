#include "catch2/catch.hpp"
#include "control/atmosphere_driver.hpp"
#include "physics/rrtmgp/rrtmgp.hpp"
#include "physics/shoc/shoc.hpp"

namespace {

TEST_CASE("scream", "stubs") {
  int val_driver  = scream::control::driver_stub();
  int val_rrtmpg  = scream::rrtmgp::rrtmgp_stub();
  int val_shoc    = scream::shoc::shoc_stub();

  REQUIRE(val_driver == 42);
  REQUIRE(val_rrtmpg == 42);
  REQUIRE(val_shoc == 42);
}

} // empty namespace
