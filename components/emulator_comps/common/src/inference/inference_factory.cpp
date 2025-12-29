/**
 * @file inference_factory.cpp
 * @brief Factory functions for creating inference backends.
 *
 * Provides centralized backend creation with compile-time conditional
 * support for different ML frameworks.
 */

#include "inference_backend.hpp"
#include "stub_backend.hpp"

#ifdef EMULATOR_HAS_LIBTORCH
#include "libtorch_backend.hpp"
#endif

namespace emulator {
namespace inference {

/**
 * @brief Create an inference backend by type.
 *
 * Creates an uninitialized backend instance. The caller must call
 * initialize() on the returned object before using it.
 *
 * @param type Backend type to create
 * @return Unique pointer to new backend, or StubBackend if type unavailable
 *
 * @note If a requested backend is not available at compile time (e.g.,
 *       LIBTORCH without EMULATOR_HAS_LIBTORCH defined), falls back to Stub.
 */
std::unique_ptr<InferenceBackend> create_backend(BackendType type) {
  switch (type) {
  case BackendType::STUB:
    return std::make_unique<StubBackend>();

#ifdef EMULATOR_HAS_LIBTORCH
  case BackendType::LIBTORCH:
    return std::make_unique<LibTorchBackend>();
#endif

  default:
    // Fall back to stub for unknown or unavailable backends
    return std::make_unique<StubBackend>();
  }
}

/**
 * @brief Create and initialize a backend from configuration.
 *
 * Convenience function that creates a backend and calls initialize().
 * If initialization fails, returns nullptr.
 *
 * @param config Configuration with backend type and parameters
 * @return Initialized backend, or nullptr on failure
 */
std::unique_ptr<InferenceBackend>
create_backend(const InferenceConfig &config) {
  auto backend = create_backend(config.backend);

  if (backend && !backend->initialize(config)) {
    return nullptr;
  }

  return backend;
}

} // namespace inference
} // namespace emulator
