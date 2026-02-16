/**
 * @file create_backend.cpp
 * @brief Factory for creating inference backends.
 */

#include "create_backend.hpp"
#include "stub_backend.hpp"

namespace emulator {
namespace inference {

std::shared_ptr<InferenceBackend>
create_backend(BackendType type, const InferenceConfig &config) {
  switch (type) {
  case BackendType::STUB:
    return std::make_shared<StubBackend>(config);
  default:
    return std::make_shared<StubBackend>(config);
  }
}

} // namespace inference
} // namespace emulator
