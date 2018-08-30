#include "catch2/catch.hpp"
#include "coupler/coupler.hpp"
#include "rrtmgp/rrtmgp.hpp"
#include "shoc/shoc.hpp"

namespace {

TEST_CASE("scream", "stubs") {
  int val_coupler = scream::coupler::coupler_stub();
  int val_rrtmpg  = scream::rrtmgp::rrtmgp_stub();
  int val_shoc    = scream::shoc::shoc_stub();

  REQUIRE(val_coupler == 42);
  REQUIRE(val_rrtmpg == 42);
  REQUIRE(val_shoc == 42);
}

}  // namespace
