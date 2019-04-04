#include "catch2/catch.hpp"
#include "physics/shoc/shoc.hpp"

namespace {

TEST_CASE("Shoc", "stub") {
  int val = scream::shoc::shoc_stub();
  REQUIRE(val == 42);
}

} // empty namespace
