/**
 * @file create_inference_backend.cpp
 * @brief Factory for creating inference backends.
 *
 * A switch rather than a registry: there are two backends, and a registry
 * that has to be populated explicitly (a static library may drop an object
 * file whose symbols nothing references, silently un-registering a backend)
 * is more machinery than two names justify.
 */

#include "create_inference_backend.hpp"

#include "inference_error.hpp"
#include "stub_inference_backend.hpp"

namespace emulator {
namespace inference {

std::vector<std::string> available_backends() {
  return {"stub"};
}

std::shared_ptr<InferenceBackend>
create_backend(const InferenceConfig &config, const InferenceContext &context) {
  std::shared_ptr<InferenceBackend> backend;

  if (config.backend.empty() || config.backend == "stub" ||
      config.backend == "none") {
    backend = std::make_shared<StubBackend>(config, context);
  } else {
    std::string names;
    for (const auto &n : available_backends()) {
      names += (names.empty() ? "" : ", ") + n;
    }
    EMULATOR_INFER_REQUIRE(false, "Unknown inference backend '"
                                      << config.backend
                                      << "'. This build has: " << names << ".");
  }

  backend->initialize();
  return backend;
}

} // namespace inference
} // namespace emulator
