/**
 * @file test_data_source.cpp
 * @brief Unit tests for DataSource interface and MctDataSource.
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "data_source.hpp"
#include "data_view.hpp"
#include "mct_data_source.hpp"

#include <vector>

namespace emulator {
namespace test {

// ── Mock DataSource ─────────────────────────────────────────────────

/**
 * @brief Trivial DataSource for testing the interface contract.
 */
class MockDataSource : public DataSource {
public:
  double fill_value = 42.0;
  bool import_called = false;
  bool export_called = false;
  std::vector<double> exported_data;

  std::vector<FieldSpec> available_fields() const override {
    return {FieldSpec("mock_a", "units_a"), FieldSpec("mock_b", "units_b")};
  }

  void import_data(DataView &view) override {
    import_called = true;
    view.fill(fill_value);
  }

  void export_data(const DataView &view) override {
    export_called = true;
    exported_data.assign(view.data(), view.data() + view.size());
  }
};

// ── DataSource interface tests ──────────────────────────────────────

TEST_CASE("DataSource mock import_data", "[data_source]") {
  MockDataSource source;
  DataView view("test", 5);
  view.add_field(FieldSpec("mock_a"));
  view.allocate();

  source.import_data(view);

  REQUIRE(source.import_called);
  REQUIRE(view(0, 0) == Approx(42.0));
  REQUIRE(view(0, 4) == Approx(42.0));
}

TEST_CASE("DataSource mock export_data", "[data_source]") {
  MockDataSource source;
  DataView view("test", 3);
  view.add_field(FieldSpec("mock_a"));
  view.add_field(FieldSpec("mock_b"));
  view.allocate();
  view.fill(7.0);

  source.export_data(view);

  REQUIRE(source.export_called);
  REQUIRE(source.exported_data.size() == 6);
  for (double v : source.exported_data) {
    REQUIRE(v == Approx(7.0));
  }
}

TEST_CASE("DataSource default export_data is no-op", "[data_source]") {
  // A base class pointer with default export_data should not throw
  class ReadOnlySource : public DataSource {
  public:
    std::vector<FieldSpec> available_fields() const override {
      return {FieldSpec("ro")};
    }
    void import_data(DataView &view) override { view.fill(1.0); }
    // export_data() uses default no-op
  };

  ReadOnlySource source;
  DataView view("test", 3);
  view.add_field(FieldSpec("ro"));
  view.allocate();
  view.fill(99.0);

  REQUIRE_NOTHROW(source.export_data(view));
}

TEST_CASE("DataSource available_fields", "[data_source]") {
  MockDataSource source;
  auto fields = source.available_fields();
  REQUIRE(fields.size() == 2);
  REQUIRE(fields[0].name == "mock_a");
  REQUIRE(fields[1].name == "mock_b");
}

// ── MctDataSource tests ─────────────────────────────────────────────

TEST_CASE("MctDataSource import FIELD_MAJOR", "[data_source][mct]") {
  const int nfields = 2;
  const int npoints = 4;

  // MCT rAttr: FIELD_MAJOR, field0={1,2,3,4}, field1={5,6,7,8}
  std::vector<double> mct = {1, 2, 3, 4, 5, 6, 7, 8};

  MctDataSource source(mct.data(), nfields, npoints, DataLayout::FIELD_MAJOR,
                       {FieldSpec("a"), FieldSpec("b")});

  DataView view("v", npoints, DataLayout::FIELD_MAJOR);
  view.add_field(FieldSpec("a"));
  view.add_field(FieldSpec("b"));
  view.allocate();

  source.import_data(view);

  REQUIRE(view(0, 0) == Approx(1.0));
  REQUIRE(view(0, 3) == Approx(4.0));
  REQUIRE(view(1, 0) == Approx(5.0));
  REQUIRE(view(1, 3) == Approx(8.0));
}

TEST_CASE("MctDataSource import with transposition", "[data_source][mct]") {
  const int nfields = 2;
  const int npoints = 3;

  // MCT: FIELD_MAJOR, field0={10,20,30}, field1={40,50,60}
  std::vector<double> mct = {10, 20, 30, 40, 50, 60};

  MctDataSource source(mct.data(), nfields, npoints, DataLayout::FIELD_MAJOR,
                       {FieldSpec("x"), FieldSpec("y")});

  // DataView is POINT_MAJOR
  DataView view("v", npoints, DataLayout::POINT_MAJOR);
  view.add_field(FieldSpec("x"));
  view.add_field(FieldSpec("y"));
  view.allocate();

  source.import_data(view);

  REQUIRE(view(0, 0) == Approx(10.0));
  REQUIRE(view(0, 2) == Approx(30.0));
  REQUIRE(view(1, 0) == Approx(40.0));
  REQUIRE(view(1, 2) == Approx(60.0));
}

TEST_CASE("MctDataSource export round-trip", "[data_source][mct]") {
  const int nfields = 2;
  const int npoints = 3;

  std::vector<double> mct_in = {1, 2, 3, 4, 5, 6};
  std::vector<double> mct_out(6, 0.0);

  MctDataSource source(mct_in.data(), nfields, npoints, DataLayout::FIELD_MAJOR,
                       {FieldSpec("a"), FieldSpec("b")});
  source.set_export_buffer(mct_out.data());

  DataView view("v", npoints, DataLayout::FIELD_MAJOR);
  view.add_field(FieldSpec("a"));
  view.add_field(FieldSpec("b"));
  view.allocate();

  source.import_data(view);
  source.export_data(view);

  // Round-trip should preserve values
  for (int i = 0; i < 6; ++i) {
    REQUIRE(mct_out[static_cast<std::size_t>(i)] ==
            Approx(mct_in[static_cast<std::size_t>(i)]));
  }
}

TEST_CASE("MctDataSource with field_map", "[data_source][mct]") {
  const int nfields = 3;
  const int npoints = 3;

  // Source: 3 fields. field0={1,2,3}, field1={4,5,6}, field2={7,8,9}
  std::vector<double> mct = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  // Only import fields 0 and 2 → dst fields 0 and 1
  std::vector<std::pair<int, int>> fmap = {{0, 0}, {2, 1}};
  MctDataSource source(mct.data(), nfields, npoints, DataLayout::FIELD_MAJOR,
                       {FieldSpec("a"), FieldSpec("c")}, fmap);

  DataView view("v", npoints, DataLayout::FIELD_MAJOR);
  view.add_field(FieldSpec("a"));
  view.add_field(FieldSpec("c"));
  view.allocate();

  source.import_data(view);

  REQUIRE(view(0, 0) == Approx(1.0));
  REQUIRE(view(0, 2) == Approx(3.0));
  REQUIRE(view(1, 0) == Approx(7.0));
  REQUIRE(view(1, 2) == Approx(9.0));
}

TEST_CASE("MctDataSource read-only (no export buffer)", "[data_source][mct]") {
  std::vector<double> mct = {1, 2, 3};
  MctDataSource source(mct.data(), 1, 3, DataLayout::FIELD_MAJOR,
                       {FieldSpec("f")});

  DataView view("v", 3);
  view.add_field(FieldSpec("f"));
  view.allocate();

  source.import_data(view);
  // export_data should be a no-op (no export buffer set)
  REQUIRE_NOTHROW(source.export_data(view));
}

TEST_CASE("MctDataSource null import buffer throws", "[data_source][mct]") {
  MctDataSource source(nullptr, 1, 3, DataLayout::FIELD_MAJOR,
                       {FieldSpec("f")});

  DataView view("v", 3);
  view.add_field(FieldSpec("f"));
  view.allocate();

  REQUIRE_THROWS_AS(source.import_data(view), std::runtime_error);
}

TEST_CASE("MctDataSource update_import_buffer", "[data_source][mct]") {
  std::vector<double> mct1 = {1, 2, 3};
  std::vector<double> mct2 = {10, 20, 30};

  MctDataSource source(mct1.data(), 1, 3, DataLayout::FIELD_MAJOR,
                       {FieldSpec("f")});

  DataView view("v", 3);
  view.add_field(FieldSpec("f"));
  view.allocate();

  source.import_data(view);
  REQUIRE(view(0, 0) == Approx(1.0));

  source.update_import_buffer(mct2.data());
  source.import_data(view);
  REQUIRE(view(0, 0) == Approx(10.0));
}

} // namespace test
} // namespace emulator
