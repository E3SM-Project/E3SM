/**
 * @file test_tensor_adapter.cpp
 * @brief Tests for TensorAdapter and RawTensorAdapter.
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "data_view.hpp"
#include "raw_tensor_adapter.hpp"

#include <cstring>
#include <numeric>
#include <vector>

namespace emulator {
namespace test {

static DataView make_view(int npoints, int nfields,
                          DataLayout layout = DataLayout::POINT_MAJOR) {
  DataView view("test", npoints, layout);
  for (int f = 0; f < nfields; ++f) {
    view.add_field(FieldSpec("f" + std::to_string(f)));
  }
  view.allocate();
  // Fill with known values
  for (int f = 0; f < nfields; ++f) {
    for (int p = 0; p < npoints; ++p) {
      view(f, p) = f * 100.0 + p;
    }
  }
  return view;
}

// ── RawTensorAdapter default shape ──────────────────────────────────

TEST_CASE("RawTensorAdapter default shape POINT_MAJOR", "[tensor]") {
  auto view = make_view(5, 3, DataLayout::POINT_MAJOR);
  RawTensorAdapter adapter;

  auto shape = adapter.tensor_shape(view);
  REQUIRE(shape.size() == 2);
  REQUIRE(shape[0] == 5); // npoints
  REQUIRE(shape[1] == 3); // nfields
}

TEST_CASE("RawTensorAdapter default shape FIELD_MAJOR", "[tensor]") {
  auto view = make_view(5, 3, DataLayout::FIELD_MAJOR);
  RawTensorAdapter adapter;

  auto shape = adapter.tensor_shape(view);
  REQUIRE(shape.size() == 2);
  REQUIRE(shape[0] == 3); // nfields
  REQUIRE(shape[1] == 5); // npoints
}

// ── RawTensorAdapter data_ptr zero-copy ─────────────────────────────

TEST_CASE("RawTensorAdapter data_ptr returns view data", "[tensor]") {
  auto view = make_view(4, 2);
  RawTensorAdapter adapter;

  const double *ptr = adapter.data_ptr(view);
  REQUIRE(ptr == view.data());
}

// ── RawTensorAdapter to_tensor / from_tensor ────────────────────────

TEST_CASE("RawTensorAdapter to_tensor copies data", "[tensor]") {
  auto view = make_view(3, 2);
  RawTensorAdapter adapter;

  std::vector<double> buffer(view.size(), 0.0);
  adapter.to_tensor(view, buffer.data());

  // Buffer should match view's raw data
  for (std::size_t i = 0; i < view.size(); ++i) {
    REQUIRE(buffer[i] == Approx(view.data()[i]));
  }
}

TEST_CASE("RawTensorAdapter from_tensor writes data", "[tensor]") {
  auto view = make_view(3, 2);
  RawTensorAdapter adapter;

  std::vector<double> buffer(view.size());
  for (std::size_t i = 0; i < buffer.size(); ++i) {
    buffer[i] = static_cast<double>(i) * 7.0;
  }

  adapter.from_tensor(buffer.data(), view);

  for (std::size_t i = 0; i < view.size(); ++i) {
    REQUIRE(view.data()[i] == Approx(buffer[i]));
  }
}

// ── RawTensorAdapter custom shape ───────────────────────────────────

TEST_CASE("RawTensorAdapter custom shape", "[tensor]") {
  // 12 elements: 4 points x 3 fields
  auto view = make_view(4, 3);
  RawTensorAdapter adapter;

  // Reshape as {1, 3, 2, 2} → batch=1, C=3, H=2, W=2
  adapter.set_custom_shape({1, 3, 2, 2});

  auto shape = adapter.tensor_shape(view);
  REQUIRE(shape == std::vector<int64_t>{1, 3, 2, 2});

  // data_ptr still works (reinterpretation only)
  REQUIRE(adapter.data_ptr(view) == view.data());
}

TEST_CASE("RawTensorAdapter custom shape size mismatch", "[tensor]") {
  auto view = make_view(4, 3); // 12 elements
  RawTensorAdapter adapter;

  adapter.set_custom_shape({1, 3, 3, 3}); // 27 != 12

  // tensor_shape should throw
  REQUIRE_THROWS_AS(adapter.tensor_shape(view), std::runtime_error);

  // data_ptr should return nullptr
  REQUIRE(adapter.data_ptr(view) == nullptr);
}

TEST_CASE("RawTensorAdapter clear custom shape", "[tensor]") {
  auto view = make_view(4, 3);
  RawTensorAdapter adapter;

  adapter.set_custom_shape({1, 3, 2, 2});
  REQUIRE(adapter.has_custom_shape());

  adapter.clear_custom_shape();
  REQUIRE_FALSE(adapter.has_custom_shape());

  auto shape = adapter.tensor_shape(view);
  REQUIRE(shape[0] == 4); // back to default npoints
  REQUIRE(shape[1] == 3);
}

// ── Error handling ──────────────────────────────────────────────────

TEST_CASE("RawTensorAdapter null buffer errors", "[tensor]") {
  auto view = make_view(2, 2);
  RawTensorAdapter adapter;

  REQUIRE_THROWS_AS(adapter.to_tensor(view, nullptr), std::runtime_error);
  REQUIRE_THROWS_AS(adapter.from_tensor(nullptr, view), std::runtime_error);
}

TEST_CASE("RawTensorAdapter empty custom shape throws", "[tensor]") {
  RawTensorAdapter adapter;
  REQUIRE_THROWS_AS(adapter.set_custom_shape({}), std::runtime_error);
}

// ── Round-trip: MCT → DataView → tensor → DataView → MCT ───────────

TEST_CASE("Full round-trip MCT-DataView-tensor", "[tensor][integration]") {
  const int nfields = 3;
  const int npoints = 5;

  // 1. Simulate MCT rAttr buffer (FIELD_MAJOR)
  std::vector<double> mct_buffer(nfields * npoints);
  for (int f = 0; f < nfields; ++f) {
    for (int p = 0; p < npoints; ++p) {
      mct_buffer[static_cast<std::size_t>(f * npoints + p)] =
          (f + 1) * 100.0 + p;
    }
  }

  // 2. Import into POINT_MAJOR DataView
  DataView view("round_trip", npoints, DataLayout::POINT_MAJOR);
  for (int f = 0; f < nfields; ++f) {
    view.add_field(FieldSpec("f" + std::to_string(f)));
  }
  view.allocate();
  view.import_from(mct_buffer.data(), nfields, npoints,
                   DataLayout::FIELD_MAJOR);

  // 3. Pass to tensor adapter (zero-copy pointer)
  RawTensorAdapter adapter;
  const double *tensor_ptr = adapter.data_ptr(view);
  REQUIRE(tensor_ptr != nullptr);

  // 4. Simulate ML: multiply all values by 2
  std::vector<double> ml_output(view.size());
  for (std::size_t i = 0; i < view.size(); ++i) {
    ml_output[i] = tensor_ptr[i] * 2.0;
  }

  // 5. Write ML output back into a different DataView (export)
  DataView export_view("export", npoints, DataLayout::POINT_MAJOR);
  for (int f = 0; f < nfields; ++f) {
    export_view.add_field(FieldSpec("f" + std::to_string(f)));
  }
  export_view.allocate();
  adapter.from_tensor(ml_output.data(), export_view);

  // 6. Export back to MCT-style FIELD_MAJOR buffer
  std::vector<double> mct_out(nfields * npoints, 0.0);
  export_view.export_to(mct_out.data(), nfields, npoints,
                        DataLayout::FIELD_MAJOR);

  // 7. Verify: every value should be 2x the original
  for (int f = 0; f < nfields; ++f) {
    for (int p = 0; p < npoints; ++p) {
      auto idx = static_cast<std::size_t>(f * npoints + p);
      double expected = ((f + 1) * 100.0 + p) * 2.0;
      REQUIRE(mct_out[idx] == Approx(expected));
    }
  }
}

} // namespace test
} // namespace emulator
