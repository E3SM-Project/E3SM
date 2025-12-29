//==============================================================================
// Test: CouplingFields
//
// Unit tests for coupling field management.
//==============================================================================

#include "coupling_fields.hpp"
#include <catch2/catch.hpp>

using namespace emulator;

TEST_CASE("CouplingFieldsBase construction", "[unit][coupling]") {
  CouplingFieldsBase fields;

  SECTION("initializes with zero counts") {
    REQUIRE(fields.num_exports == 0);
    REQUIRE(fields.num_imports == 0);
  }
}

TEST_CASE("CouplingFieldsBase field parsing", "[unit][coupling]") {
  CouplingFieldsBase fields;

  SECTION("parses colon-separated export field list") {
    fields.initialize("a2x_Sa_z:a2x_Sa_tbot:a2x_Sa_pbot", "");

    REQUIRE(fields.num_exports == 3);
    REQUIRE(fields.get_export_index("a2x_Sa_z") == 0);
    REQUIRE(fields.get_export_index("a2x_Sa_tbot") == 1);
    REQUIRE(fields.get_export_index("a2x_Sa_pbot") == 2);
  }

  SECTION("parses colon-separated import field list") {
    fields.initialize("", "x2a_Sx_t:x2a_Faxx_sen:x2a_Faxx_lat");

    REQUIRE(fields.num_imports == 3);
    REQUIRE(fields.get_import_index("x2a_Sx_t") == 0);
    REQUIRE(fields.get_import_index("x2a_Faxx_sen") == 1);
    REQUIRE(fields.get_import_index("x2a_Faxx_lat") == 2);
  }

  SECTION("assigns sequential indices") {
    fields.initialize("field_a:field_b:field_c", "");

    REQUIRE(fields.get_export_index("field_a") == 0);
    REQUIRE(fields.get_export_index("field_b") == 1);
    REQUIRE(fields.get_export_index("field_c") == 2);
  }
}

TEST_CASE("CouplingFieldsBase field lookup", "[unit][coupling]") {
  CouplingFieldsBase fields;

  // Initialize with some fields
  fields.initialize("a2x_Sa_z:a2x_Faxa_rainc", "x2a_Sx_t:x2a_Faxx_sen");

  SECTION("finds registered export fields") {
    REQUIRE(fields.get_export_index("a2x_Sa_z") >= 0);
    REQUIRE(fields.get_export_index("a2x_Faxa_rainc") >= 0);
  }

  SECTION("finds registered import fields") {
    REQUIRE(fields.get_import_index("x2a_Sx_t") >= 0);
    REQUIRE(fields.get_import_index("x2a_Faxx_sen") >= 0);
  }

  SECTION("returns -1 for unknown fields") {
    REQUIRE(fields.get_export_index("nonexistent") == -1);
    REQUIRE(fields.get_import_index("nonexistent") == -1);
  }
}

TEST_CASE("CouplingFieldsBase edge cases", "[unit][coupling]") {
  CouplingFieldsBase fields;

  SECTION("handles empty export string") {
    fields.initialize("", "x2a_Sx_t");
    REQUIRE(fields.num_exports == 0);
    REQUIRE(fields.num_imports == 1);
  }

  SECTION("handles empty import string") {
    fields.initialize("a2x_Sa_z", "");
    REQUIRE(fields.num_exports == 1);
    REQUIRE(fields.num_imports == 0);
  }

  SECTION("handles both empty") {
    fields.initialize("", "");
    REQUIRE(fields.num_exports == 0);
    REQUIRE(fields.num_imports == 0);
  }

  SECTION("handles single field") {
    fields.initialize("a2x_Sa_z", "");
    REQUIRE(fields.num_exports == 1);
    REQUIRE(fields.get_export_index("a2x_Sa_z") == 0);
  }
}
