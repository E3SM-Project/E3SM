#include <catch2/catch.hpp>

#include "ekat/scream_pack.hpp"

extern "C"{
int test_array_io ();
}

namespace {

TEST_CASE("array_io_mod", "test_array_io") {

  int nerr = test_array_io();

  REQUIRE (nerr==0);
}

} // empty namespace
