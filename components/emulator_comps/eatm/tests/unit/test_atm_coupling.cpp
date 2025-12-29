//==============================================================================
// Test: AtmCoupling
//
// Unit tests for atmosphere coupling utilities.
//==============================================================================

#include "impl/atm_coupling.hpp"
#include "impl/atm_field_manager.hpp"
#include <catch2/catch.hpp>

using namespace emulator;
using namespace emulator::impl;

TEST_CASE("AtmCouplingIndices default values", "[unit][eatm][coupling]") {
  AtmCouplingIndices indices;

  SECTION("export indices default to -1") {
    REQUIRE(indices.Sa_z == -1);
    REQUIRE(indices.Sa_u == -1);
    REQUIRE(indices.Sa_tbot == -1);
    REQUIRE(indices.Faxa_rainc == -1);
  }

  SECTION("import indices default to -1") {
    REQUIRE(indices.Sx_t == -1);
    REQUIRE(indices.Faxx_sen == -1);
    REQUIRE(indices.Faxx_lat == -1);
    REQUIRE(indices.Sf_ifrac == -1);
  }
}

TEST_CASE("AtmCouplingIndices initialization", "[unit][eatm][coupling]") {
  AtmCouplingIndices indices;
  CouplingFieldsBase fields;

  SECTION("initializes from CouplingFieldsBase") {
    // Setup coupling fields
    fields.initialize("a2x_Sa_z:a2x_Sa_u:a2x_Sa_v:a2x_Sa_tbot", // exports
                      "x2a_Sx_t:x2a_Faxx_sen:x2a_Faxx_lat"      // imports
    );

    indices.initialize(fields);

    // After initialization, indices should be set based on field positions
    // The actual behavior depends on how initialize() maps field names
  }
}

TEST_CASE("CouplingFieldsBase parsing", "[unit][eatm][coupling]") {
  CouplingFieldsBase fields;

  SECTION("parses colon-separated export list") {
    fields.initialize("a2x_Sa_z:a2x_Sa_u:a2x_Sa_tbot", "");

    REQUIRE(fields.num_exports == 3);
    REQUIRE(fields.get_export_index("a2x_Sa_z") == 0);
    REQUIRE(fields.get_export_index("a2x_Sa_u") == 1);
    REQUIRE(fields.get_export_index("a2x_Sa_tbot") == 2);
  }

  SECTION("parses colon-separated import list") {
    fields.initialize("", "x2a_Sx_t:x2a_Faxx_sen");

    REQUIRE(fields.num_imports == 2);
    REQUIRE(fields.get_import_index("x2a_Sx_t") == 0);
    REQUIRE(fields.get_import_index("x2a_Faxx_sen") == 1);
  }

  SECTION("returns -1 for unknown fields") {
    fields.initialize("a2x_Sa_z", "x2a_Sx_t");

    REQUIRE(fields.get_export_index("unknown") == -1);
    REQUIRE(fields.get_import_index("unknown") == -1);
  }

  SECTION("handles empty string") {
    fields.initialize("", "");

    REQUIRE(fields.num_exports == 0);
    REQUIRE(fields.num_imports == 0);
  }
}

TEST_CASE("import_atm_fields function", "[unit][eatm][coupling]") {
  const int ncol = 10;
  const int nfields = 5;

  AtmFieldManager field_mgr;
  field_mgr.allocate(ncol);

  AtmCouplingIndices indices;

  // Create mock import buffer
  std::vector<double> import_data(nfields * ncol, 0.0);

  SECTION("handles zero indices gracefully") {
    // With all indices at -1, import should be a no-op
    import_atm_fields(import_data.data(), ncol, nfields, indices, field_mgr);
    SUCCEED("No crash with uninitialized indices");
  }
}

TEST_CASE("export_atm_fields function", "[unit][eatm][coupling]") {
  const int ncol = 10;
  const int nfields = 5;

  AtmFieldManager field_mgr;
  field_mgr.allocate(ncol);
  field_mgr.set_defaults(ncol);

  AtmCouplingIndices indices;

  // Create export buffer
  std::vector<double> export_data(nfields * ncol, 0.0);

  SECTION("handles zero indices gracefully") {
    // With all indices at -1, export should be a no-op
    export_atm_fields(export_data.data(), ncol, nfields, indices, field_mgr);
    SUCCEED("No crash with uninitialized indices");
  }
}
