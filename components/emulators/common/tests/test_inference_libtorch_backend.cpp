// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "create_inference_backend.hpp"
#include "inference_error.hpp"

#include <numeric>
#include <string>
#include <vector>

namespace emulator {
namespace inference {
namespace test {

namespace {

using Dims = std::vector<std::int64_t>;

// The shape the fixtures were traced on: [batch, channels, ny, nx]
const Dims kDims{1, 3, 2, 4};

// Fixtures are written at build time by fixtures/make_libtorch_fixture.py
InferenceConfig fixture(const std::string &name) {
  InferenceConfig config;
  config.model_path = std::string(EMULATOR_TEST_TORCH_FIXTURE_DIR) + "/" + name;
  return config;
}

std::shared_ptr<InferenceBackend> create(const InferenceConfig &config) {
  return create_backend(BackendType::LIBTORCH, config);
}

// 1, 2, ..., n
std::vector<double> ramp(std::size_t n = 24) {
  std::vector<double> v(n);
  std::iota(v.begin(), v.end(), 1.0);
  return v;
}

bool infer(InferenceBackend &backend, const std::vector<double> &x,
           const Dims &x_dims, std::vector<double> &y, const Dims &y_dims) {
  TensorMap inputs;
  inputs.wrap("x", x.data(), x_dims);
  TensorMap outputs;
  outputs.wrap("y", y.data(), y_dims);
  return backend.infer(inputs, outputs);
}

} // namespace

TEST_CASE("LibTorchBackend runs a TorchScript module", "[libtorch]") {
  auto backend = create(fixture("affine.pt"));
  REQUIRE(backend->name() == "LibTorch");

  const auto x = ramp();
  std::vector<double> y(24, -999.0);

  // Exact despite the double -> float32 -> double round trip
  REQUIRE(infer(*backend, x, kDims, y, kDims));
  for (std::size_t i = 0; i < y.size(); ++i) {
    REQUIRE(y[i] == 2.0 * x[i] + 1.0);
  }

  // Repeated calls work, and the input view is never written
  y.assign(24, -999.0);
  REQUIRE(infer(*backend, x, kDims, y, kDims));
  REQUIRE(y[23] == 49.0);
  REQUIRE(x == ramp());
}

TEST_CASE("LibTorchBackend flat-array inference", "[libtorch]") {
  auto config = fixture("affine.pt");
  config.input_channels = 3;
  config.output_channels = 3;
  auto backend = create(config);

  const auto x = ramp(6);
  std::vector<double> y(6, 0.0);

  REQUIRE(backend->infer(x.data(), y.data(), 2));
  REQUIRE(y[5] == 13.0);
}

TEST_CASE("LibTorchBackend dtype sets the precision used", "[libtorch]") {
  // 0.1 is not representable, so float32 and float64 results differ
  const double exact = 2.0 * 0.1 + 1.0;
  const std::vector<double> x{0.1};
  std::vector<double> y{0.0};
  const Dims one{1, 1, 1, 1};

  SECTION("float32 is the default") {
    auto backend = create(fixture("affine.pt"));
    REQUIRE(infer(*backend, x, one, y, one));
    REQUIRE(y[0] == Approx(exact));
    REQUIRE(y[0] != exact);
  }

  SECTION("float64 is exact") {
    auto config = fixture("affine.pt");
    config.set("dtype", "float64");
    auto backend = create(config);
    REQUIRE(infer(*backend, x, one, y, one));
    REQUIRE(y[0] == exact);
  }
}

TEST_CASE("LibTorchBackend fills several outputs in order", "[libtorch]") {
  const std::string name = GENERATE("tuple.pt", "list.pt");
  auto backend = create(fixture(name));

  const auto x = ramp();
  std::vector<double> a(24, 0.0);
  std::vector<double> b(24, 0.0);

  TensorMap inputs;
  inputs.wrap("x", x.data(), kDims);
  TensorMap outputs;
  outputs.wrap("a", a.data(), kDims);
  outputs.wrap("b", b.data(), kDims);

  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(a[23] == 2.0 * x[23] + 1.0);
  REQUIRE(b[23] == x[23] - 1.0);

  SECTION("and refuses a different number of outputs") {
    TensorMap one_output;
    one_output.wrap("a", a.data(), kDims);
    REQUIRE_THROWS_WITH(backend->infer(inputs, one_output),
                        Catch::Contains("2 tensor(s) but 1 output(s)"));
  }
}

TEST_CASE("LibTorchBackend output shape may differ from input", "[libtorch]") {
  auto backend = create(fixture("channels.pt"));

  // Channel 0 is 1..8, channel 1 is 9..16, channel 2 is 17..24
  const auto x = ramp();
  std::vector<double> y(8, 0.0);

  REQUIRE(infer(*backend, x, kDims, y, {1, 1, 2, 4}));
  for (std::size_t i = 0; i < 8; ++i) {
    REQUIRE(y[i] == x[i] + 10.0 * x[i + 8] + 100.0 * x[i + 16]);
  }
}

TEST_CASE("LibTorchBackend reports shape errors", "[libtorch]") {
  SECTION("an output that does not match what the module returned") {
    auto backend = create(fixture("affine.pt"));
    const auto x = ramp();
    std::vector<double> y(24, 0.0);

    // Same element count, different shape: never reshaped to fit
    REQUIRE_THROWS_WITH(infer(*backend, x, kDims, y, {1, 3, 4, 2}),
                        Catch::Contains("[1, 3, 2, 4]") &&
                            Catch::Contains("y[1,3,4,2]"));
  }

  SECTION("an input the module rejects") {
    auto backend = create(fixture("channels.pt"));
    const auto x = ramp(32);
    std::vector<double> y(8, 0.0);

    REQUIRE_THROWS_WITH(infer(*backend, x, {1, 4, 2, 4}, y, {1, 1, 2, 4}),
                        Catch::Contains("forward()") &&
                            Catch::Contains("x[1,4,2,4]"));
  }

  SECTION("no inputs at all") {
    auto backend = create(fixture("affine.pt"));
    TensorMap inputs;
    TensorMap outputs;
    REQUIRE_THROWS_AS(backend->infer(inputs, outputs), InferenceError);
  }
}

TEST_CASE("LibTorchBackend accepts empty tensors", "[libtorch]") {
  // e.g. an MPI rank that owns no columns
  auto backend = create(fixture("affine.pt"));

  const std::vector<double> x;
  std::vector<double> y;
  REQUIRE(infer(*backend, x, {0, 3}, y, {0, 3}));
}

TEST_CASE("LibTorchBackend reports a model it cannot load", "[libtorch]") {
  SECTION("no path") {
    REQUIRE_THROWS_WITH(create(InferenceConfig()),
                        Catch::Contains("model_path"));
  }

  SECTION("a missing file") {
    REQUIRE_THROWS_WITH(create(fixture("no_such_module.pt")),
                        Catch::Contains("no_such_module.pt"));
  }

  SECTION("a file that is not TorchScript") {
    InferenceConfig config;
    config.model_path = __FILE__;
    REQUIRE_THROWS_AS(create(config), InferenceError);
  }
}

TEST_CASE("LibTorchBackend never swaps CUDA for the CPU", "[libtorch]") {
  auto config = fixture("affine.pt");
  config.set("device", "cuda");

  // Either there is a GPU and it is used, or creation fails
  std::shared_ptr<InferenceBackend> backend;
  try {
    backend = create(config);
  } catch (const InferenceError &e) {
    REQUIRE_THAT(e.what(), Catch::Contains("no CUDA device"));
    return;
  }

  const auto x = ramp();
  std::vector<double> y(24, 0.0);
  REQUIRE(infer(*backend, x, kDims, y, kDims));
  REQUIRE(y[23] == 2.0 * x[23] + 1.0);
}

TEST_CASE("LibTorchBackend validates its options", "[libtorch]") {
  auto config = fixture("affine.pt");

  SECTION("device") {
    config.set("device", "tpu");
    REQUIRE_THROWS_WITH(create(config), Catch::Contains("cpu, cuda or cuda:N"));

    config.set("device", "cuda:99");
    REQUIRE_THROWS_AS(create(config), InferenceError);
  }

  SECTION("dtype") {
    config.set("dtype", "bfloat16");
    REQUIRE_THROWS_WITH(create(config), Catch::Contains("float32 or float64"));
  }

  SECTION("num_threads") {
    config.set("num_threads", "-1");
    REQUIRE_THROWS_AS(create(config), InferenceError);

    config.set("num_threads", "2");
    auto backend = create(config);
    const auto x = ramp();
    std::vector<double> y(24, 0.0);
    REQUIRE(infer(*backend, x, kDims, y, kDims));
    REQUIRE(y[0] == 3.0);
  }
}

TEST_CASE("LibTorchBackend refuses inference after finalize", "[libtorch]") {
  auto backend = create(fixture("affine.pt"));
  backend->finalize();
  backend->finalize();

  const auto x = ramp();
  std::vector<double> y(24, 0.0);
  REQUIRE_THROWS_WITH(infer(*backend, x, kDims, y, kDims),
                      Catch::Contains("finalized"));
}

} // namespace test
} // namespace inference
} // namespace emulator
