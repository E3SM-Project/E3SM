/**
 * @file test_io_adapter.cpp
 * @brief Tests for the IOAdapter abstract interface using a mock.
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "data_view.hpp"
#include "decomp_info.hpp"
#include "io_adapter.hpp"

#include <string>
#include <vector>

namespace emulator {
namespace test {

/**
 * @brief Mock IOAdapter that records calls and metadata for verification.
 *
 * Does not perform actual IO — just validates that the correct
 * DataView metadata is passed through the interface.
 */
class MockIOAdapter : public IOAdapter {
public:
  // Recorded state
  bool init_called = false;
  bool define_called = false;
  bool write_called = false;
  bool read_called = false;
  bool close_called = false;

  int recorded_nfields = 0;
  int recorded_npoints = 0;
  int recorded_nglobal = 0;
  std::vector<std::string> recorded_field_names;
  std::vector<int64_t> recorded_dof;
  std::vector<int> recorded_global_dims;
  std::string recorded_filename;

  void init_io(const DataView &view) override {
    init_called = true;
    recorded_nfields = view.num_fields();
    recorded_npoints = view.num_points();
    if (view.has_decomp()) {
      recorded_nglobal = view.decomp().npoints_global;
      recorded_dof = view.decomp().dof;
      recorded_global_dims = view.decomp().global_dims;
    }
  }

  void define_variables(const DataView &view) override {
    define_called = true;
    recorded_field_names.clear();
    for (int i = 0; i < view.num_fields(); ++i) {
      recorded_field_names.push_back(view.field_spec(i).name);
    }
  }

  void write(const DataView &view, const std::string &filename) override {
    write_called = true;
    recorded_filename = filename;
    // Verify we can access field data
    for (int i = 0; i < view.num_fields(); ++i) {
      (void)view.field_data(i); // should not throw
    }
  }

  void read(DataView &view, const std::string &filename) override {
    read_called = true;
    recorded_filename = filename;
    // Simulate reading: fill with known values
    for (int f = 0; f < view.num_fields(); ++f) {
      for (int p = 0; p < view.num_points(); ++p) {
        view(f, p) = static_cast<double>(f * 1000 + p);
      }
    }
  }

  void close() override { close_called = true; }
};

// ── init_io passes decomposition metadata ───────────────────────────

TEST_CASE("IOAdapter init_io receives decomp metadata", "[io_adapter]") {
  DataView view("io_test", 10, DataLayout::FIELD_MAJOR);
  view.add_field(FieldSpec("Sa_z", "m", "Bottom height", 100));
  view.add_field(FieldSpec("Sa_u", "m/s", "Zonal wind", 100));
  view.allocate();

  // Attach decomp for IO
  std::vector<int64_t> dof(10);
  for (int i = 0; i < 10; ++i)
    dof[static_cast<std::size_t>(i)] = i + 1;
  view.set_decomp(DecompInfo(10, 100, dof, {100}));

  MockIOAdapter io;
  io.init_io(view);

  REQUIRE(io.init_called);
  REQUIRE(io.recorded_nfields == 2);
  REQUIRE(io.recorded_npoints == 10);
  REQUIRE(io.recorded_nglobal == 100);
  REQUIRE(io.recorded_dof.size() == 10);
  REQUIRE(io.recorded_dof[0] == 1);
  REQUIRE(io.recorded_dof[9] == 10);
  REQUIRE(io.recorded_global_dims == std::vector<int>{100});
}

// ── init_io works without decomp ────────────────────────────────────

TEST_CASE("IOAdapter init_io without decomp", "[io_adapter]") {
  DataView view("io_test", 5);
  view.add_field(FieldSpec("f"));
  view.allocate();

  MockIOAdapter io;
  io.init_io(view);

  REQUIRE(io.init_called);
  REQUIRE(io.recorded_nfields == 1);
  REQUIRE(io.recorded_npoints == 5);
  REQUIRE(io.recorded_nglobal == 0); // no decomp
}

// ── define_variables receives field specs ────────────────────────────

TEST_CASE("IOAdapter define_variables receives field names", "[io_adapter]") {
  DataView view("io_test", 5);
  view.add_field(FieldSpec("Sa_z"));
  view.add_field(FieldSpec("Sa_u"));
  view.add_field(FieldSpec("Sa_v"));
  view.allocate();

  MockIOAdapter io;
  io.define_variables(view);

  REQUIRE(io.define_called);
  REQUIRE(io.recorded_field_names.size() == 3);
  REQUIRE(io.recorded_field_names[0] == "Sa_z");
  REQUIRE(io.recorded_field_names[1] == "Sa_u");
  REQUIRE(io.recorded_field_names[2] == "Sa_v");
}

// ── write passes data correctly ─────────────────────────────────────

TEST_CASE("IOAdapter write receives filename and data", "[io_adapter]") {
  DataView view("io_test", 4, DataLayout::FIELD_MAJOR);
  view.add_field(FieldSpec("f0"));
  view.add_field(FieldSpec("f1"));
  view.allocate();
  view.fill(99.0);

  MockIOAdapter io;
  io.write(view, "/tmp/output.nc");

  REQUIRE(io.write_called);
  REQUIRE(io.recorded_filename == "/tmp/output.nc");
}

// ── read modifies DataView ──────────────────────────────────────────

TEST_CASE("IOAdapter read populates DataView", "[io_adapter]") {
  DataView view("io_test", 3);
  view.add_field(FieldSpec("a"));
  view.add_field(FieldSpec("b"));
  view.allocate();
  view.zero();

  MockIOAdapter io;
  io.read(view, "/tmp/input.nc");

  REQUIRE(io.read_called);
  REQUIRE(io.recorded_filename == "/tmp/input.nc");

  // Mock fills with f*1000+p
  REQUIRE(view(0, 0) == Approx(0.0));
  REQUIRE(view(0, 2) == Approx(2.0));
  REQUIRE(view(1, 0) == Approx(1000.0));
  REQUIRE(view(1, 2) == Approx(1002.0));
}

// ── close is callable ───────────────────────────────────────────────

TEST_CASE("IOAdapter close", "[io_adapter]") {
  MockIOAdapter io;
  io.close();
  REQUIRE(io.close_called);
}

// ── Full lifecycle ──────────────────────────────────────────────────

TEST_CASE("IOAdapter full lifecycle", "[io_adapter][integration]") {
  DataView view("lifecycle", 5, DataLayout::FIELD_MAJOR);
  view.add_field(FieldSpec("Sa_z", "m", "Bottom height", 50));
  view.add_field(FieldSpec("Sa_tbot", "K", "Temperature", 50));
  view.allocate();
  view.fill(273.15);

  MockIOAdapter io;

  // 1. Init IO
  io.init_io(view);
  REQUIRE(io.init_called);

  // 2. Define variables
  io.define_variables(view);
  REQUIRE(io.define_called);
  REQUIRE(io.recorded_field_names[0] == "Sa_z");

  // 3. Write
  io.write(view, "/scratch/run/output.nc");
  REQUIRE(io.write_called);

  // 4. Read into a different view
  DataView view2("read_target", 5, DataLayout::FIELD_MAJOR);
  view2.add_field(FieldSpec("Sa_z"));
  view2.add_field(FieldSpec("Sa_tbot"));
  view2.allocate();

  io.read(view2, "/scratch/run/input.nc");
  REQUIRE(io.read_called);

  // 5. Close
  io.close();
  REQUIRE(io.close_called);
}

} // namespace test
} // namespace emulator
