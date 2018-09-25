#include "catch2/catch.hpp"

#include "share/scream_config.hpp"
#include "share/scream_pack.hpp"

extern "C"{
int test_array_io ();
}

namespace {

TEST_CASE("array_io_mod", "test_array_io") {

  int nerr = test_array_io();

  REQUIRE (nerr==0);
}

} // empty namespace

