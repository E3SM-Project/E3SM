/**
 * @file test_data_pipeline.cpp
 * @brief Integration tests for DataPipeline at increasing coupling levels.
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "data_pipeline.hpp"
#include "data_view.hpp"
#include "io_adapter.hpp"
#include "mct_data_source.hpp"
#include "raw_tensor_adapter.hpp"

#include <cstring>
#include <vector>

namespace emulator {
namespace test {

// ── Re-usable MockIOAdapter ─────────────────────────────────────────

class MockIOAdapter : public IOAdapter {
public:
  bool write_called = false;
  bool read_called = false;
  std::string recorded_filename;
  double read_fill_value = 99.0;

  void init_io(const DataView &) override {}
  void define_variables(const DataView &) override {}
  void write(const DataView &, const std::string &f) override {
    write_called = true;
    recorded_filename = f;
  }
  void read(DataView &view, const std::string &f) override {
    read_called = true;
    recorded_filename = f;
    view.fill(read_fill_value);
  }
  void close() override {}
};

// ── Pipeline construction ───────────────────────────────────────────

TEST_CASE("DataPipeline basic construction", "[pipeline]") {
  DataPipeline pipe("test", 10);

  REQUIRE(pipe.name() == "test");
  REQUIRE_FALSE(pipe.is_initialized());
}

TEST_CASE("DataPipeline initialize with import fields only", "[pipeline]") {
  DataPipeline pipe("import_only", 5);
  pipe.add_import_field(FieldSpec("a"));
  pipe.add_import_field(FieldSpec("b"));
  pipe.initialize();

  REQUIRE(pipe.is_initialized());
  REQUIRE(pipe.import_view().num_fields() == 2);
  REQUIRE(pipe.import_view().num_points() == 5);
  REQUIRE(pipe.import_view().name() == "import_only_import");
}

TEST_CASE("DataPipeline initialize with export fields only", "[pipeline]") {
  DataPipeline pipe("export_only", 5);
  pipe.add_export_field(FieldSpec("x"));
  pipe.initialize();

  REQUIRE(pipe.is_initialized());
  REQUIRE(pipe.export_view().num_fields() == 1);
  REQUIRE(pipe.export_view().name() == "export_only_export");
}

TEST_CASE("DataPipeline initialize with no fields throws", "[pipeline]") {
  DataPipeline pipe("empty", 5);
  REQUIRE_THROWS_AS(pipe.initialize(), std::runtime_error);
}

TEST_CASE("DataPipeline double initialize throws", "[pipeline]") {
  DataPipeline pipe("test", 5);
  pipe.add_import_field(FieldSpec("f"));
  pipe.initialize();
  REQUIRE_THROWS_AS(pipe.initialize(), std::runtime_error);
}

// ── Pipeline + MctDataSource ────────────────────────────────────────

TEST_CASE("DataPipeline ingest from MctDataSource",
          "[pipeline][mct][integration]") {
  const int nfields = 2;
  const int npoints = 4;

  // MCT: FIELD_MAJOR, field0={1,2,3,4}, field1={5,6,7,8}
  std::vector<double> mct_in = {1, 2, 3, 4, 5, 6, 7, 8};

  DataPipeline pipe("atm", npoints);
  pipe.set_import_source(std::make_unique<MctDataSource>(
      mct_in.data(), nfields, npoints, DataLayout::FIELD_MAJOR,
      std::vector<FieldSpec>{FieldSpec("a"), FieldSpec("b")}));
  pipe.add_import_field(FieldSpec("a"));
  pipe.add_import_field(FieldSpec("b"));
  pipe.initialize();

  pipe.ingest();

  auto &view = pipe.import_view();
  REQUIRE(view(0, 0) == Approx(1.0));
  REQUIRE(view(0, 3) == Approx(4.0));
  REQUIRE(view(1, 0) == Approx(5.0));
  REQUIRE(view(1, 3) == Approx(8.0));
}

TEST_CASE("DataPipeline egest to MctDataSource",
          "[pipeline][mct][integration]") {
  const int nfields = 2;
  const int npoints = 3;

  std::vector<double> mct_out(6, 0.0);

  auto source = std::make_unique<MctDataSource>(
      nullptr, nfields, npoints, DataLayout::FIELD_MAJOR,
      std::vector<FieldSpec>{FieldSpec("x"), FieldSpec("y")});
  source->set_export_buffer(mct_out.data());

  DataPipeline pipe("atm", npoints);
  pipe.set_export_source(std::move(source));
  pipe.add_export_field(FieldSpec("x"));
  pipe.add_export_field(FieldSpec("y"));
  // Need at least one view to initialize
  pipe.add_import_field(FieldSpec("dummy"));
  pipe.initialize();

  // Write known values into export view
  auto &ev = pipe.export_view();
  ev(0, 0) = 10.0;
  ev(0, 1) = 20.0;
  ev(0, 2) = 30.0;
  ev(1, 0) = 40.0;
  ev(1, 1) = 50.0;
  ev(1, 2) = 60.0;

  pipe.egest();

  // Verify MCT buffer received the correct values
  // DataView is POINT_MAJOR, MCT expects FIELD_MAJOR
  // field0={10,20,30}, field1={40,50,60}
  REQUIRE(mct_out[0] == Approx(10.0));
  REQUIRE(mct_out[1] == Approx(20.0));
  REQUIRE(mct_out[2] == Approx(30.0));
  REQUIRE(mct_out[3] == Approx(40.0));
  REQUIRE(mct_out[4] == Approx(50.0));
  REQUIRE(mct_out[5] == Approx(60.0));
}

// ── Pipeline + RawTensorAdapter ─────────────────────────────────────

TEST_CASE("DataPipeline tensor access", "[pipeline][tensor]") {
  const int npoints = 4;
  const int nfields = 3;

  DataPipeline pipe("ml", npoints);
  pipe.set_tensor_adapter(std::make_unique<RawTensorAdapter>());
  for (int f = 0; f < nfields; ++f) {
    pipe.add_import_field(FieldSpec("f" + std::to_string(f)));
    pipe.add_export_field(FieldSpec("f" + std::to_string(f)));
  }
  pipe.initialize();

  // Fill import view
  auto &iv = pipe.import_view();
  for (int f = 0; f < nfields; ++f) {
    for (int p = 0; p < npoints; ++p) {
      iv(f, p) = f * 10.0 + p;
    }
  }

  // Read tensor shape
  auto shape = pipe.tensor_shape();
  REQUIRE(shape.size() == 2);
  REQUIRE(shape[0] == npoints);
  REQUIRE(shape[1] == nfields);

  // Zero-copy tensor access
  const double *tdata = pipe.tensor_data();
  REQUIRE(tdata == iv.data());

  // Simulate ML: multiply by 2, write to export view
  std::vector<double> ml_out(iv.size());
  for (std::size_t i = 0; i < iv.size(); ++i) {
    ml_out[i] = tdata[i] * 2.0;
  }
  pipe.write_tensor(ml_out.data());

  // Verify export view
  auto &ev = pipe.export_view();
  REQUIRE(ev(0, 0) == Approx(0.0));  // 0 * 2
  REQUIRE(ev(0, 3) == Approx(6.0));  // 3 * 2
  REQUIRE(ev(2, 0) == Approx(40.0)); // 20 * 2
}

// ── Pipeline + MockIOAdapter ────────────────────────────────────────

TEST_CASE("DataPipeline IO write and read", "[pipeline][io]") {
  DataPipeline pipe("io_test", 5);
  auto io = std::make_unique<MockIOAdapter>();
  auto *io_ptr = io.get();

  pipe.set_io_adapter(std::move(io));
  pipe.add_import_field(FieldSpec("f0"));
  pipe.initialize();

  pipe.import_view().fill(42.0);
  pipe.write_io("/tmp/test.nc");

  REQUIRE(io_ptr->write_called);
  REQUIRE(io_ptr->recorded_filename == "/tmp/test.nc");

  pipe.read_io("/tmp/input.nc");
  REQUIRE(io_ptr->read_called);
  REQUIRE(io_ptr->recorded_filename == "/tmp/input.nc");
  // Mock fills with 99.0
  REQUIRE(pipe.import_view()(0, 0) == Approx(99.0));
}

// ── Pipeline error handling ─────────────────────────────────────────

TEST_CASE("DataPipeline missing components throw", "[pipeline]") {
  DataPipeline pipe("err", 5);
  pipe.add_import_field(FieldSpec("f"));
  pipe.initialize();

  REQUIRE_THROWS_AS(pipe.ingest(), std::runtime_error);
  REQUIRE_THROWS_AS(pipe.egest(), std::runtime_error);
  REQUIRE_THROWS_AS(pipe.tensor_data(), std::runtime_error);
  REQUIRE_THROWS_AS(pipe.tensor_shape(), std::runtime_error);
  REQUIRE_THROWS_AS(pipe.write_tensor(nullptr), std::runtime_error);
  REQUIRE_THROWS_AS(pipe.write_io("f"), std::runtime_error);
  REQUIRE_THROWS_AS(pipe.read_io("f"), std::runtime_error);
}

// ── Full integration: MCT → DataView → tensor → infer → export ─────

TEST_CASE("Full pipeline: MCT to ML round-trip",
          "[pipeline][mct][tensor][integration]") {
  const int nfields = 3;
  const int npoints = 5;

  // 1. Simulate MCT import buffer (FIELD_MAJOR)
  std::vector<double> mct_in(nfields * npoints);
  for (int f = 0; f < nfields; ++f) {
    for (int p = 0; p < npoints; ++p) {
      mct_in[static_cast<std::size_t>(f * npoints + p)] = (f + 1) * 100.0 + p;
    }
  }

  // 2. Prepare MCT export buffer
  std::vector<double> mct_out(nfields * npoints, 0.0);

  // 3. Build pipeline
  std::vector<FieldSpec> fields;
  for (int f = 0; f < nfields; ++f) {
    fields.push_back(FieldSpec("f" + std::to_string(f)));
  }

  auto import_src = std::make_unique<MctDataSource>(
      mct_in.data(), nfields, npoints, DataLayout::FIELD_MAJOR, fields);

  auto export_src = std::make_unique<MctDataSource>(
      nullptr, nfields, npoints, DataLayout::FIELD_MAJOR, fields);
  export_src->set_export_buffer(mct_out.data());

  DataPipeline pipe("atm", npoints);
  pipe.set_import_source(std::move(import_src));
  pipe.set_export_source(std::move(export_src));
  pipe.set_tensor_adapter(std::make_unique<RawTensorAdapter>());

  for (const auto &f : fields) {
    pipe.add_import_field(f);
    pipe.add_export_field(f);
  }
  pipe.initialize();

  // 4. Ingest MCT data
  pipe.ingest();

  // 5. Get tensor, simulate ML (multiply by 2)
  const double *tensor = pipe.tensor_data();
  std::vector<double> ml_result(pipe.import_view().size());
  for (std::size_t i = 0; i < ml_result.size(); ++i) {
    ml_result[i] = tensor[i] * 2.0;
  }

  // 6. Write ML output to export view
  pipe.write_tensor(ml_result.data());

  // 7. Egest to MCT
  pipe.egest();

  // 8. Verify: MCT output should be 2x the input
  for (int f = 0; f < nfields; ++f) {
    for (int p = 0; p < npoints; ++p) {
      auto idx = static_cast<std::size_t>(f * npoints + p);
      double expected = ((f + 1) * 100.0 + p) * 2.0;
      REQUIRE(mct_out[idx] == Approx(expected));
    }
  }
}

// ── Pipeline-free composable test ───────────────────────────────────

TEST_CASE("Composable: DataView + adapter without pipeline",
          "[composable][integration]") {
  const int npoints = 4;
  const int nfields = 2;

  // Direct DataView use — no pipeline needed
  DataView view("standalone", npoints, DataLayout::POINT_MAJOR);
  view.add_field(FieldSpec("temp"));
  view.add_field(FieldSpec("pressure"));
  view.allocate();

  // Fill manually
  for (int f = 0; f < nfields; ++f) {
    for (int p = 0; p < npoints; ++p) {
      view(f, p) = (f + 1) * 10.0 + p;
    }
  }

  // Use tensor adapter directly
  RawTensorAdapter adapter;
  const double *ptr = adapter.data_ptr(view);
  REQUIRE(ptr == view.data());

  auto shape = adapter.tensor_shape(view);
  REQUIRE(shape[0] == npoints);
  REQUIRE(shape[1] == nfields);
}

// ── auto_register_fields ────────────────────────────────────────────

TEST_CASE("DataPipeline auto_register_fields", "[pipeline]") {
  const int npoints = 3;
  std::vector<double> mct = {1, 2, 3, 4, 5, 6};

  auto src = std::make_unique<MctDataSource>(
      mct.data(), 2, npoints, DataLayout::FIELD_MAJOR,
      std::vector<FieldSpec>{FieldSpec("auto_a"), FieldSpec("auto_b")});

  DataPipeline pipe("auto", npoints);
  pipe.set_import_source(std::move(src));
  pipe.auto_register_fields();
  pipe.initialize();

  REQUIRE(pipe.import_view().num_fields() == 2);
  REQUIRE(pipe.import_view().field_spec(0).name == "auto_a");
  REQUIRE(pipe.import_view().field_spec(1).name == "auto_b");

  pipe.ingest();
  REQUIRE(pipe.import_view()(0, 0) == Approx(1.0));
  REQUIRE(pipe.import_view()(1, 0) == Approx(4.0));
}

} // namespace test
} // namespace emulator
