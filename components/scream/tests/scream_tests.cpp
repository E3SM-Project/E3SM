#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "coupler/coupler.hpp"
#include "p3/p3.hpp"
#include "rrtmgp/rrtmgp.hpp"
#include "shoc/shoc.hpp"

namespace {

TEST_CASE("Scream", "stubs") {
  int val_coupler = Scream::Coupler::coupler_stub();
  int val_p3      = Scream::P3::p3_stub();
  int val_rrtmpg  = Scream::Rrtmgp::rrtmgp_stub();
  int val_shoc    = Scream::Shoc::shoc_stub();

  REQUIRE(val_coupler == 42);
  REQUIRE(val_p3 == 42);
  REQUIRE(val_rrtmpg == 42);
  REQUIRE(val_shoc == 42);
}

} // empty namespace
