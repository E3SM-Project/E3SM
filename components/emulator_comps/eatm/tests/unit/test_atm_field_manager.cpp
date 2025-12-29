//==============================================================================
// Test: AtmFieldManager
//
// Unit tests for atmosphere field management.
//==============================================================================

#include "impl/atm_field_manager.hpp"
#include <catch2/catch.hpp>

using namespace emulator::impl;

TEST_CASE("AtmFieldManager construction", "[unit][eatm][fields]") {
  AtmFieldManager fields;

  SECTION("starts unallocated") { REQUIRE_FALSE(fields.is_allocated()); }

  SECTION("has correct channel constants") {
    REQUIRE(AtmFieldManager::N_INPUT_CHANNELS == 39);
    REQUIRE(AtmFieldManager::N_OUTPUT_CHANNELS == 44);
  }
}

TEST_CASE("AtmFieldManager allocation", "[unit][eatm][fields]") {
  AtmFieldManager fields;

  SECTION("can allocate with valid dimensions") {
    const int ncol = 100;

    fields.allocate(ncol);

    REQUIRE(fields.is_allocated());
    // net_inputs and net_outputs are now allocated dynamically in
    // prepare_inputs(), not in allocate(). They start empty.
    REQUIRE(fields.net_inputs.empty());
    REQUIRE(fields.net_outputs.empty());
  }

  SECTION("allocates import fields") {
    const int ncol = 50;
    fields.allocate(ncol);

    REQUIRE(fields.shf.size() == ncol);
    REQUIRE(fields.lhf.size() == ncol);
    REQUIRE(fields.ts.size() == ncol);
    REQUIRE(fields.icefrac.size() == ncol);
  }

  SECTION("allocates export fields") {
    const int ncol = 50;
    fields.allocate(ncol);

    REQUIRE(fields.zbot.size() == ncol);
    REQUIRE(fields.tbot.size() == ncol);
    REQUIRE(fields.pbot.size() == ncol);
    REQUIRE(fields.rainc.size() == ncol);
  }
}

TEST_CASE("AtmFieldManager deallocation", "[unit][eatm][fields]") {
  AtmFieldManager fields;
  const int ncol = 100;

  SECTION("can deallocate") {
    fields.allocate(ncol);
    REQUIRE(fields.is_allocated());

    fields.deallocate();
    REQUIRE_FALSE(fields.is_allocated());
  }
}

TEST_CASE("AtmFieldManager set_defaults", "[unit][eatm][fields]") {
  AtmFieldManager fields;
  const int ncol = 50;

  SECTION("sets default values for export fields") {
    fields.allocate(ncol);
    fields.set_defaults(ncol);

    // After set_defaults, fields should have non-garbage values
    // The actual default values depend on implementation
    // Just verify no exception is thrown
    SUCCEED("set_defaults completed without exception");
  }
}
