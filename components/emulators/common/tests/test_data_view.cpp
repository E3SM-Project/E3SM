/**
 * @file test_data_view.cpp
 * @brief Tests for DataView, FieldSpec, and DecompInfo.
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "data_view.hpp"
#include "decomp_info.hpp"
#include "field_spec.hpp"

#include <cmath>
#include <vector>

namespace emulator {
namespace test {

// ── FieldSpec ───────────────────────────────────────────────────────

TEST_CASE("FieldSpec default construction", "[field_spec]") {
  FieldSpec fs;
  REQUIRE(fs.name.empty());
  REQUIRE(fs.units.empty());
  REQUIRE(fs.long_name.empty());
  REQUIRE(fs.global_size == 0);
}

TEST_CASE("FieldSpec parameterized construction", "[field_spec]") {
  FieldSpec fs("Sa_z", "m", "Bottom height", 100);
  REQUIRE(fs.name == "Sa_z");
  REQUIRE(fs.units == "m");
  REQUIRE(fs.long_name == "Bottom height");
  REQUIRE(fs.global_size == 100);
}

// ── DecompInfo ──────────────────────────────────────────────────────

TEST_CASE("DecompInfo validation", "[decomp_info]") {
  SECTION("valid: dof matches npoints_local") {
    DecompInfo d(10, 100, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    REQUIRE_NOTHROW(d.validate());
  }

  SECTION("valid: empty dof is ok") {
    DecompInfo d(10, 100);
    REQUIRE_NOTHROW(d.validate());
  }

  SECTION("invalid: dof size mismatch") {
    DecompInfo d(10, 100, {1, 2, 3}); // 3 != 10
    REQUIRE_THROWS_AS(d.validate(), std::runtime_error);
  }
}

// ── DataView construction ───────────────────────────────────────────

TEST_CASE("DataView construction", "[data_view]") {
  DataView view("test", 5, DataLayout::POINT_MAJOR);

  REQUIRE(view.name() == "test");
  REQUIRE(view.num_points() == 5);
  REQUIRE(view.num_fields() == 0);
  REQUIRE_FALSE(view.is_allocated());
  REQUIRE(view.layout() == DataLayout::POINT_MAJOR);
  REQUIRE_FALSE(view.has_decomp());
}

TEST_CASE("DataView with optional DecompInfo", "[data_view]") {
  DataView view("decomp_test", 10);

  REQUIRE_FALSE(view.has_decomp());
  REQUIRE_THROWS_AS(view.decomp(), std::runtime_error);

  // Attach decomp
  std::vector<int64_t> dof(10);
  for (int i = 0; i < 10; ++i)
    dof[static_cast<std::size_t>(i)] = i + 1;
  view.set_decomp(DecompInfo(10, 100, dof, {100}));

  REQUIRE(view.has_decomp());
  REQUIRE(view.decomp().npoints_global == 100);
  REQUIRE(view.decomp().dof.size() == 10);
}

TEST_CASE("DataView decomp mismatch throws", "[data_view]") {
  DataView view("mismatch", 5);
  REQUIRE_THROWS_AS(view.set_decomp(DecompInfo(10, 100)), std::runtime_error);
}

// ── Field registration ──────────────────────────────────────────────

TEST_CASE("DataView field registration", "[data_view]") {
  DataView view("v", 4);

  view.add_field(FieldSpec("Sa_z", "m"));
  view.add_field(FieldSpec("Sa_u", "m/s"));
  view.add_field(FieldSpec("Sa_v", "m/s"));

  REQUIRE(view.num_fields() == 3);
  REQUIRE(view.field_index("Sa_z") == 0);
  REQUIRE(view.field_index("Sa_u") == 1);
  REQUIRE(view.field_index("Sa_v") == 2);
  REQUIRE(view.field_spec("Sa_z").units == "m");

  SECTION("duplicate name throws") {
    REQUIRE_THROWS_AS(view.add_field(FieldSpec("Sa_z")), std::runtime_error);
  }

  SECTION("empty name throws") {
    REQUIRE_THROWS_AS(view.add_field(FieldSpec("")), std::runtime_error);
  }
}

TEST_CASE("DataView add_fields batch", "[data_view]") {
  DataView view("v", 4);

  std::vector<FieldSpec> specs = {FieldSpec("a"), FieldSpec("b"),
                                  FieldSpec("c")};
  view.add_fields(specs);

  REQUIRE(view.num_fields() == 3);
}

// ── Allocation ──────────────────────────────────────────────────────

TEST_CASE("DataView allocation", "[data_view]") {
  DataView view("v", 3);

  view.add_field(FieldSpec("f1"));
  view.add_field(FieldSpec("f2"));

  REQUIRE_FALSE(view.is_allocated());

  view.allocate();

  REQUIRE(view.is_allocated());
  REQUIRE(view.size() == 6); // 2 fields * 3 points
  REQUIRE(view.data() != nullptr);

  SECTION("double allocate throws") {
    REQUIRE_THROWS_AS(view.allocate(), std::runtime_error);
  }

  SECTION("add_field after allocate throws") {
    REQUIRE_THROWS_AS(view.add_field(FieldSpec("f3")), std::runtime_error);
  }
}

TEST_CASE("DataView allocate with no fields throws", "[data_view]") {
  DataView view("v", 3);
  REQUIRE_THROWS_AS(view.allocate(), std::runtime_error);
}

// ── Element access (POINT_MAJOR) ────────────────────────────────────

TEST_CASE("DataView element access POINT_MAJOR", "[data_view]") {
  DataView view("v", 3, DataLayout::POINT_MAJOR);

  view.add_field(FieldSpec("f0"));
  view.add_field(FieldSpec("f1"));
  view.allocate();

  // Write via operator()
  for (int f = 0; f < 2; ++f) {
    for (int p = 0; p < 3; ++p) {
      view(f, p) = f * 10.0 + p;
    }
  }

  // Read back
  REQUIRE(view(0, 0) == Approx(0.0));
  REQUIRE(view(0, 2) == Approx(2.0));
  REQUIRE(view(1, 0) == Approx(10.0));
  REQUIRE(view(1, 2) == Approx(12.0));

  // Raw buffer: POINT_MAJOR → data[point * nfields + field]
  const double *raw = view.data();
  REQUIRE(raw[0 * 2 + 0] == Approx(0.0));  // p=0, f=0
  REQUIRE(raw[0 * 2 + 1] == Approx(10.0)); // p=0, f=1
  REQUIRE(raw[2 * 2 + 0] == Approx(2.0));  // p=2, f=0
  REQUIRE(raw[2 * 2 + 1] == Approx(12.0)); // p=2, f=1

  // shape
  auto [d0, d1] = view.shape();
  REQUIRE(d0 == 3); // npoints
  REQUIRE(d1 == 2); // nfields
}

// ── Element access (FIELD_MAJOR) ────────────────────────────────────

TEST_CASE("DataView element access FIELD_MAJOR", "[data_view]") {
  DataView view("v", 3, DataLayout::FIELD_MAJOR);

  view.add_field(FieldSpec("f0"));
  view.add_field(FieldSpec("f1"));
  view.allocate();

  for (int f = 0; f < 2; ++f) {
    for (int p = 0; p < 3; ++p) {
      view(f, p) = f * 10.0 + p;
    }
  }

  // Raw buffer: FIELD_MAJOR → data[field * npoints + point]
  const double *raw = view.data();
  REQUIRE(raw[0 * 3 + 0] == Approx(0.0));  // f=0, p=0
  REQUIRE(raw[0 * 3 + 2] == Approx(2.0));  // f=0, p=2
  REQUIRE(raw[1 * 3 + 0] == Approx(10.0)); // f=1, p=0
  REQUIRE(raw[1 * 3 + 2] == Approx(12.0)); // f=1, p=2

  // field_data for FIELD_MAJOR returns contiguous slice
  const double *fd0 = view.field_data(0);
  REQUIRE(fd0[0] == Approx(0.0));
  REQUIRE(fd0[1] == Approx(1.0));
  REQUIRE(fd0[2] == Approx(2.0));

  const double *fd1 = view.field_data("f1");
  REQUIRE(fd1[0] == Approx(10.0));
  REQUIRE(fd1[1] == Approx(11.0));
  REQUIRE(fd1[2] == Approx(12.0));

  auto [d0, d1] = view.shape();
  REQUIRE(d0 == 2); // nfields
  REQUIRE(d1 == 3); // npoints
}

// ── Import/Export same layout ───────────────────────────────────────

TEST_CASE("DataView import_from same layout", "[data_view]") {
  DataView view("v", 4, DataLayout::FIELD_MAJOR);

  view.add_field(FieldSpec("a"));
  view.add_field(FieldSpec("b"));
  view.allocate();

  // Simulate MCT rAttr: FIELD_MAJOR, 2 fields, 4 points
  // field 0: {1, 2, 3, 4}, field 1: {5, 6, 7, 8}
  std::vector<double> src = {1, 2, 3, 4, 5, 6, 7, 8};

  view.import_from(src.data(), 2, 4, DataLayout::FIELD_MAJOR);

  REQUIRE(view(0, 0) == Approx(1.0));
  REQUIRE(view(0, 3) == Approx(4.0));
  REQUIRE(view(1, 0) == Approx(5.0));
  REQUIRE(view(1, 3) == Approx(8.0));
}

// ── Import/Export with layout transposition ─────────────────────────

TEST_CASE("DataView import_from with transposition", "[data_view]") {
  // DataView is POINT_MAJOR, source is FIELD_MAJOR (like MCT rAttr)
  DataView view("v", 3, DataLayout::POINT_MAJOR);
  view.add_field(FieldSpec("x"));
  view.add_field(FieldSpec("y"));
  view.allocate();

  // MCT-style: FIELD_MAJOR, 2 fields, 3 points
  // field 0: {10, 20, 30}, field 1: {40, 50, 60}
  std::vector<double> src = {10, 20, 30, 40, 50, 60};

  view.import_from(src.data(), 2, 3, DataLayout::FIELD_MAJOR);

  REQUIRE(view(0, 0) == Approx(10.0));
  REQUIRE(view(0, 1) == Approx(20.0));
  REQUIRE(view(0, 2) == Approx(30.0));
  REQUIRE(view(1, 0) == Approx(40.0));
  REQUIRE(view(1, 1) == Approx(50.0));
  REQUIRE(view(1, 2) == Approx(60.0));

  // Export back to FIELD_MAJOR buffer
  std::vector<double> dst(6, 0.0);
  view.export_to(dst.data(), 2, 3, DataLayout::FIELD_MAJOR);

  REQUIRE(dst[0] == Approx(10.0));
  REQUIRE(dst[1] == Approx(20.0));
  REQUIRE(dst[2] == Approx(30.0));
  REQUIRE(dst[3] == Approx(40.0));
  REQUIRE(dst[4] == Approx(50.0));
  REQUIRE(dst[5] == Approx(60.0));
}

// ── Import/Export with field_map ─────────────────────────────────────

TEST_CASE("DataView import_from with field_map", "[data_view]") {
  DataView view("v", 3, DataLayout::FIELD_MAJOR);

  view.add_field(FieldSpec("selected_a"));
  view.add_field(FieldSpec("selected_c"));
  view.allocate();

  // Source has 3 fields (a, b, c), we only want fields 0 and 2
  // FIELD_MAJOR: field0={1,2,3}, field1={4,5,6}, field2={7,8,9}
  std::vector<double> src = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::vector<std::pair<int, int>> fmap = {{0, 0}, {2, 1}};
  view.import_from(src.data(), 3, 3, DataLayout::FIELD_MAJOR, fmap);

  REQUIRE(view(0, 0) == Approx(1.0)); // src field 0 → dst field 0
  REQUIRE(view(0, 2) == Approx(3.0));
  REQUIRE(view(1, 0) == Approx(7.0)); // src field 2 → dst field 1
  REQUIRE(view(1, 2) == Approx(9.0));
}

// ── Zero and fill ───────────────────────────────────────────────────

TEST_CASE("DataView zero and fill", "[data_view]") {
  DataView view("v", 4);
  view.add_field(FieldSpec("f"));
  view.allocate();

  view.fill(42.0);
  for (std::size_t i = 0; i < view.size(); ++i) {
    REQUIRE(view.data()[i] == Approx(42.0));
  }

  view.zero();
  for (std::size_t i = 0; i < view.size(); ++i) {
    REQUIRE(view.data()[i] == Approx(0.0));
  }
}

// ── copy_from ───────────────────────────────────────────────────────

TEST_CASE("DataView copy_from same layout", "[data_view]") {
  DataView a("a", 3, DataLayout::FIELD_MAJOR);
  DataView b("b", 3, DataLayout::FIELD_MAJOR);

  a.add_field(FieldSpec("f0"));
  a.add_field(FieldSpec("f1"));
  a.allocate();

  b.add_field(FieldSpec("g0"));
  b.add_field(FieldSpec("g1"));
  b.allocate();

  a(0, 0) = 1.0;
  a(0, 1) = 2.0;
  a(1, 0) = 3.0;
  a(1, 1) = 4.0;

  b.copy_from(a);

  REQUIRE(b(0, 0) == Approx(1.0));
  REQUIRE(b(0, 1) == Approx(2.0));
  REQUIRE(b(1, 0) == Approx(3.0));
  REQUIRE(b(1, 1) == Approx(4.0));
}

TEST_CASE("DataView copy_from transposing layout", "[data_view]") {
  DataView a("a", 3, DataLayout::FIELD_MAJOR);
  DataView b("b", 3, DataLayout::POINT_MAJOR);

  a.add_field(FieldSpec("f0"));
  a.add_field(FieldSpec("f1"));
  a.allocate();

  b.add_field(FieldSpec("g0"));
  b.add_field(FieldSpec("g1"));
  b.allocate();

  a(0, 0) = 10.0;
  a(0, 2) = 30.0;
  a(1, 1) = 50.0;

  b.copy_from(a);

  // Values should be the same via (field, point) semantics
  REQUIRE(b(0, 0) == Approx(10.0));
  REQUIRE(b(0, 2) == Approx(30.0));
  REQUIRE(b(1, 1) == Approx(50.0));
}

// ── Error conditions ────────────────────────────────────────────────

TEST_CASE("DataView access before allocate throws", "[data_view]") {
  DataView view("v", 3);
  view.add_field(FieldSpec("f"));

  REQUIRE_THROWS_AS(view.data(), std::runtime_error);
  REQUIRE_THROWS_AS(view.field_data(0), std::runtime_error);
  REQUIRE_THROWS_AS(view.import_from(nullptr, 1, 3, DataLayout::FIELD_MAJOR),
                    std::runtime_error);
}

TEST_CASE("DataView field_index unknown name throws", "[data_view]") {
  DataView view("v", 3);
  view.add_field(FieldSpec("f"));
  REQUIRE_THROWS_AS(view.field_index("nonexistent"), std::runtime_error);
}

} // namespace test
} // namespace emulator
