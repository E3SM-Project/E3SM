// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "inference_config.hpp"
#include "inference_error.hpp"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace emulator {
namespace inference {
namespace test {

namespace {

/// Write a temporary file and remove it when the test ends.
class ScopedFile {
public:
  ScopedFile(std::string path, const std::string &contents)
      : m_path(std::move(path)) {
    std::ofstream ofs(m_path);
    ofs << contents;
  }
  ~ScopedFile() { std::remove(m_path.c_str()); }
  const std::string &path() const { return m_path; }

private:
  std::string m_path;
};

} // namespace

TEST_CASE("InferenceConfig defaults", "[inference_config]") {
  InferenceConfig config;

  REQUIRE(config.backend == "stub");
  REQUIRE(config.input_channels == 0);
  REQUIRE(config.output_channels == 0);
  REQUIRE_FALSE(config.verbose);
  REQUIRE(config.options.empty());
}

TEST_CASE("apply() sets the named fields", "[inference_config]") {
  InferenceConfig config;
  config.apply("backend", " python ");
  config.apply("model_path", "'/scratch/ace.tar'");
  config.apply("input", "T");
  config.apply("input", "q");
  config.apply("output", "dT");
  config.apply("input_channels", "73");
  config.apply("verbose", ".true.");

  REQUIRE(config.backend == "python");
  REQUIRE(config.model_path == "/scratch/ace.tar"); // quotes stripped
  REQUIRE(config.inputs == std::vector<std::string>{"T", "q"});
  REQUIRE(config.outputs == std::vector<std::string>{"dT"});
  REQUIRE(config.input_channels == 73);
  REQUIRE(config.verbose);
}

TEST_CASE("A malformed value is refused rather than guessed at",
          "[inference_config]") {
  InferenceConfig config;
  REQUIRE_THROWS_AS(config.apply("verbose", "maybe"), InferenceError);
  REQUIRE_THROWS_AS(config.apply("input_channels", "lots"), InferenceError);
}

TEST_CASE("Unknown keys become options rather than errors",
          "[inference_config]") {
  // The Python side owns its own settings; adding one must not need a C++
  // change.
  InferenceConfig config;
  config.apply("ace_h", "4");
  config.apply("option.ace_w", "8");

  REQUIRE(config.get("ace_h") == "4");
  REQUIRE(config.get_int("ace_w", 0) == 8); // the `option.` prefix is optional
  REQUIRE(config.get("absent", "fallback") == "fallback");
  REQUIRE(config.get_bool("absent", true));
}

TEST_CASE("from_file reads a component namelist", "[inference_config]") {
  ScopedFile file("emulator_test_atm_in",
                  "nx: 90\n"
                  "ny: 45\n"
                  "! a comment line\n"
                  "inference.backend: python       # trailing comment\n"
                  "inference.emulator: ace\n"
                  "inference.model_path: /scratch/ace.tar\n"
                  "inference.input: air_temperature_0\n"
                  "inference.output: air_temperature_0\n"
                  "inference.ace_h: 4\n"
                  "not_a_setting\n");

  const auto config = InferenceConfig::from_file(file.path(), "inference.");

  REQUIRE(config.backend == "python");
  REQUIRE(config.model_path == "/scratch/ace.tar");
  REQUIRE(config.get("emulator") == "ace");
  REQUIRE(config.get_int("ace_h", 0) == 4);
  REQUIRE(config.inputs.size() == 1);
  REQUIRE(config.outputs.size() == 1);
  // Keys outside the prefix belong to the component, not to us.
  REQUIRE(config.get("nx", "unset") == "unset");
}

TEST_CASE("A missing file leaves the defaults", "[inference_config]") {
  const auto config = InferenceConfig::from_file("no_such_file_at_all");
  REQUIRE(config.backend == "stub");
}

} // namespace test
} // namespace inference
} // namespace emulator
