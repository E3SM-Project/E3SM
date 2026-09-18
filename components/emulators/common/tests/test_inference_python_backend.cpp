// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "create_inference_backend.hpp"
#include "inference_error.hpp"

#include <string>
#include <vector>

namespace emulator {
namespace inference {
namespace test {

namespace {

// An embedded interpreter does not see a virtual environment's packages,
// so the build passes the location of the ones it found
const std::string kPythonPath =
    std::string(EMULATOR_TEST_FIXTURE_DIR) + ":" + EMULATOR_TEST_SITE_PACKAGES;

// The model in fixtures/python_fixture.py
InferenceConfig fixture() {
  InferenceConfig config;
  config.model_path = "fixture.ckpt";
  config.set("python_module", "python_fixture");
  config.set("python_path", kPythonPath);
  config.set("scale", "3.0");
  return config;
}

std::shared_ptr<InferenceBackend> create(const InferenceConfig &config) {
  return create_backend(BackendType::PYTHON, config);
}

} // namespace

TEST_CASE("PythonBackend runs a Python model in place", "[python]") {
  auto backend = create(fixture());
  REQUIRE(backend->name() == "Python");

  // A writable view is still handed to the model read-only as an input;
  // the fixture raises if it is not
  std::vector<double> x{1, 2, 3, 4, 5, 6};
  std::vector<double> y(6, 0.0);
  std::vector<double> row_sum(2, 0.0);

  TensorMap inputs;
  inputs.wrap("x", x.data(), {2, 3});
  TensorMap outputs;
  outputs.wrap("y", y.data(), {2, 3});
  outputs.wrap("row_sum", row_sum.data(), {2});

  // y = scale * x + step
  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(y == std::vector<double>{4, 7, 10, 13, 16, 19});
  // Row-major: the last dimension is contiguous
  REQUIRE(row_sum == std::vector<double>{6, 15});

  // The model object persists between calls
  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(y[0] == 5.0);
  REQUIRE(x == std::vector<double>{1, 2, 3, 4, 5, 6});

  backend->finalize();
}

#ifdef EMULATOR_TEST_HAVE_TORCH
TEST_CASE("PythonBackend runs a torch.nn module", "[python][torch]") {
  // fixtures/python_torch_fixture.py: y = ReLU(W x + b), on a GPU if present
  InferenceConfig config;
  config.set("python_module", "python_torch_fixture");
  config.set("python_path", kPythonPath);
  config.verbose = true;
  auto backend = create(config);

  const std::vector<double> x{1, 2, 3, -1, -2, -3};
  std::vector<double> y(4, -999.0);

  TensorMap inputs;
  inputs.wrap("x", x.data(), {2, 3});
  TensorMap outputs;
  outputs.wrap("y", y.data(), {2, 2});

  // W = [[1,2,3],[4,5,6]], b = [0.5,-0.5]; the second row is clipped by ReLU
  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(y == std::vector<double>{14.5, 31.5, 0.0, 0.0});

  backend->finalize();
}
#endif

TEST_CASE("PythonBackend accepts empty tensors", "[python]") {
  // e.g. an MPI rank that owns no columns
  auto backend = create(fixture());

  const std::vector<double> x;
  std::vector<double> y;
  std::vector<double> row_sum;

  TensorMap inputs;
  inputs.wrap("x", x.data(), {0, 3});
  TensorMap outputs;
  outputs.wrap("y", y.data(), {0, 3});
  outputs.wrap("row_sum", row_sum.data(), {0});

  REQUIRE(backend->infer(inputs, outputs));
}

TEST_CASE("PythonBackend refuses a read-only output", "[python]") {
  auto backend = create(fixture());

  const std::vector<double> x(3, 1.0);
  TensorMap inputs;
  inputs.wrap("x", x.data(), {1, 3});
  TensorMap outputs;
  outputs.wrap("y", x.data(), {1, 3});

  REQUIRE_THROWS_WITH(backend->infer(inputs, outputs),
                      Catch::Contains("read-only"));
}

TEST_CASE("PythonBackend reports Python errors with a traceback", "[python]") {
  auto config = fixture();

  SECTION("no module named") {
    config.options.erase("python_module");
    REQUIRE_THROWS_WITH(create(config), Catch::Contains("python_module"));
  }

  SECTION("a module that cannot be imported") {
    config.set("python_module", "no_such_module");
    REQUIRE_THROWS_WITH(create(config), Catch::Contains("no_such_module") &&
                                            Catch::Contains("python_path"));
  }

  SECTION("a missing factory") {
    config.set("python_factory", "no_such_factory");
    REQUIRE_THROWS_WITH(create(config), Catch::Contains("no_such_factory"));
  }

  SECTION("a factory that raises") {
    config.set("python_factory", "create_broken");
    REQUIRE_THROWS_WITH(create(config),
                        Catch::Contains("deliberate failure in the factory") &&
                            Catch::Contains("python_fixture.py"));
  }

  SECTION("a model that raises") {
    config.set("scale", "-1");
    auto backend = create(config);

    const std::vector<double> x(3, 1.0);
    std::vector<double> y(3, 0.0);
    TensorMap inputs;
    inputs.wrap("x", x.data(), {1, 3});
    TensorMap outputs;
    outputs.wrap("y", y.data(), {1, 3});

    REQUIRE_THROWS_WITH(backend->infer(inputs, outputs),
                        Catch::Contains("deliberate failure in infer"));
  }
}

TEST_CASE("PythonBackend refuses inference after finalize", "[python]") {
  auto backend = create(fixture());
  // The fixture raises if its finalize() is called twice
  backend->finalize();
  backend->finalize();

  TensorMap inputs;
  TensorMap outputs;
  REQUIRE_THROWS_WITH(backend->infer(inputs, outputs),
                      Catch::Contains("finalized"));

  // The interpreter outlives the backend, so another one can be created
  REQUIRE(create(fixture()) != nullptr);
}

} // namespace test
} // namespace inference
} // namespace emulator
