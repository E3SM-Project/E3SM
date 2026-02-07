/**
 * @file inference_factory.cpp
 * @brief Factory for creating inference backends.
 */

#include "inference_factory.hpp"
#include "stub_backend.hpp"

namespace emulator {
namespace inference {

std::unique_ptr<InferenceBackend> create_backend(BackendType type) {
  switch (type) {
  case BackendType::STUB:
    return std::make_unique<StubBackend>();
  default:
    return std::make_unique<StubBackend>();
  }
}

} // namespace inference
} // namespace emulator
