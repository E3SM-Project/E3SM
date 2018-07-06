#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "coupler.hpp"

namespace {

TEST_CASE("Coupler", "stub") {
  int val = Scream::Coupler::coupler_stub();
  REQUIRE(val == 42);
}

} // empty namespace
