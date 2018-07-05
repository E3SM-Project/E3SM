#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "p3.hpp"

namespace {

TEST_CASE("P3", "stub") {
  int val = Scream::p3_stub();
  REQUIRE(val == 42);
}

} // empty namespace
