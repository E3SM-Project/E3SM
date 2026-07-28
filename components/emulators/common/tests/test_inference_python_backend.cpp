// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "create_inference_backend.hpp"
#include "inference_error.hpp"
#include "python_interpreter.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace emulator {
namespace inference {
namespace test {

namespace {

/// Settings pointing at tests/fixtures/emulator_fixture.py.
InferenceConfig fixture_config(const std::string &report) {
  InferenceConfig config;
  config.backend = "python";
  config.set("python_module", "emulator_fixture");
  config.set("python_path", EMULATOR_TEST_FIXTURE_DIR);
  config.set("report_path", report);
  config.set("scale", "3.0");
  config.model_path = "/not/read/by/the/fixture";
  return config;
}

std::string read_file(const std::string &path) {
  std::ifstream ifs(path);
  std::ostringstream oss;
  oss << ifs.rdbuf();
  return oss.str();
}

} // namespace

TEST_CASE("The python backend runs a model in this process", "[python]") {
  auto backend = create_backend(fixture_config("emulator_fixture_step.txt"),
                                InferenceContext());
  REQUIRE(backend->name() == "Python");

  std::vector<double> x{1.0, 2.0, 3.0, 4.0};
  std::vector<double> y(4, 0.0);

  TensorMap inputs;
  inputs.wrap("x", static_cast<const double *>(x.data()), {4});
  TensorMap outputs;
  outputs.wrap("y", y.data(), {4});

  // The fixture computes scale * x + step, writing through a view of `y`.
  // It also asserts that the input arrived read-only and the output
  // writeable, so a successful call is itself part of the assertion.
  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(y[0] == Approx(3.0 * 1.0 + 1));
  REQUIRE(y[3] == Approx(3.0 * 4.0 + 1));

  // Step two proves the emulator object persists between calls, which is
  // what an autoregressive model relies on.
  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(y[0] == Approx(3.0 * 1.0 + 2));

  // And nothing wrote back through the input view.
  REQUIRE(x[0] == 1.0);

  backend->finalize();
  std::remove("emulator_fixture_step.txt");
}

TEST_CASE("The factory receives the config and the context", "[python]") {
  const int gids[3] = {1, 5, 9};
  const double lat[3] = {-45.0, 0.0, 45.0};
  const double lon[3] = {0.0, 120.0, 240.0};

  InferenceContext context;
  context.set_grid(8, 6, 48, gids, lat, lon, 3);

  auto backend =
      create_backend(fixture_config("emulator_fixture_context.txt"), context);

  const std::string report = read_file("emulator_fixture_context.txt");
  REQUIRE(report.find("scale=3.0") != std::string::npos);
  REQUIRE(report.find("model_path=/not/read/by/the/fixture") !=
          std::string::npos);
  REQUIRE(report.find("nx=8 ny=6") != std::string::npos);
  REQUIRE(report.find("gids=[1, 5, 9]") != std::string::npos);

  backend->finalize();
  std::remove("emulator_fixture_context.txt");
}

TEST_CASE("An empty tensor still reaches the model", "[python]") {
  // A rank owning no columns is normal on a large layout, and infer() is
  // collective, so it must survive a zero-length field rather than throwing.
  auto backend = create_backend(fixture_config("emulator_fixture_empty.txt"),
                                InferenceContext());

  TensorMap inputs;
  inputs.wrap("x", static_cast<const double *>(nullptr), {0});
  TensorMap outputs;
  outputs.wrap("y", static_cast<double *>(nullptr), {0});
  REQUIRE(backend->infer(inputs, outputs));

  backend->finalize();
  std::remove("emulator_fixture_empty.txt");
}

TEST_CASE("A missing module is reported, not guessed at", "[python]") {
  InferenceConfig config;
  config.backend = "python";
  config.set("python_module", "no_such_emulator_module");
  REQUIRE_THROWS_WITH(create_backend(config, InferenceContext()),
                      Catch::Contains("python_path"));
}

TEST_CASE("A module with no factory is reported", "[python]") {
  auto config = fixture_config("emulator_fixture_nofactory.txt");
  config.set("python_factory", "not_a_function");
  REQUIRE_THROWS_AS(create_backend(config, InferenceContext()), InferenceError);
}

TEST_CASE("An error inside the model surfaces with its traceback", "[python]") {
  auto config = fixture_config("emulator_fixture_broken.txt");
  config.set("python_module", "emulator_broken");
  REQUIRE_THROWS_WITH(create_backend(config, InferenceContext()),
                      Catch::Contains("deliberate") &&
                          Catch::Contains("emulator_broken.py"));
}

TEST_CASE("A failed initialization rolls back what it acquired", "[python]") {
  // initialize() marks the backend initialized only once init_impl() returns,
  // so a throw part way through leaves finalize() declining to run: the
  // interpreter customer leaks, and the Python references taken so far are
  // released later by the destructor with no GIL held -- undefined behaviour
  // in a process where somebody else owns the interpreter.  Every failure
  // mode has to unwind to exactly where it started.
  auto &interpreter = PyInterpreter::instance();
  const int before = interpreter.num_customers();

  SECTION("the module does not import") {
    auto config = fixture_config("emulator_fixture_rollback.txt");
    config.set("python_module", "no_such_emulator_module");
    REQUIRE_THROWS(create_backend(config, InferenceContext()));
  }
  SECTION("the factory raises") {
    auto config = fixture_config("emulator_fixture_rollback.txt");
    config.set("python_module", "emulator_broken");
    REQUIRE_THROWS(create_backend(config, InferenceContext()));
  }
  SECTION("the factory is not callable") {
    auto config = fixture_config("emulator_fixture_rollback.txt");
    config.set("python_factory", "not_a_function");
    REQUIRE_THROWS(create_backend(config, InferenceContext()));
  }

  REQUIRE(interpreter.num_customers() == before);

  // And the interpreter is still usable afterwards, which is the thing a
  // leaked customer or a botched teardown would eventually cost.
  auto backend =
      create_backend(fixture_config("emulator_fixture_after.txt"),
                     InferenceContext());
  std::vector<double> x{2.0};
  std::vector<double> y{0.0};
  TensorMap inputs;
  inputs.wrap("x", static_cast<const double *>(x.data()), {1});
  TensorMap outputs;
  outputs.wrap("y", y.data(), {1});
  REQUIRE(backend->infer(inputs, outputs));
  REQUIRE(y[0] == Approx(3.0 * 2.0 + 1));
  backend->finalize();
  REQUIRE(interpreter.num_customers() == before);
  std::remove("emulator_fixture_after.txt");
  std::remove("emulator_fixture_rollback.txt");
}

} // namespace test
} // namespace inference
} // namespace emulator
