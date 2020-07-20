#include "catch2/catch.hpp"
#include "physics/rrtmgp/rrtmgp.hpp"

namespace {

TEST_CASE("rrtmgp", "stub") {
  int val = scream::rrtmgp::rrtmgp_stub();
  REQUIRE(val == 42);
}

} // empty namespace
