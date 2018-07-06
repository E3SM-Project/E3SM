#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "shoc.hpp"

namespace {

TEST_CASE("Shoc", "stub") {
  int val = Scream::SHOC::shoc_stub();
  REQUIRE(val == 42);
}

} // empty namespace
