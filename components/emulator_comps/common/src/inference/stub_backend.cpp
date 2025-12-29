/**
 * @file stub_backend.cpp
 * @brief Stub inference backend implementation.
 *
 * Provides a no-op backend for testing the emulator framework without
 * requiring ML dependencies. All outputs are set to zero.
 */

#include "stub_backend.hpp"
#include <cstring>

namespace emulator {
namespace inference {

/**
 * @brief Initialize the stub backend.
 *
 * Stores configuration for later use. Always succeeds.
 *
 * @param config Configuration (only input/output channels are used)
 * @return Always returns true
 */
bool StubBackend::initialize(const InferenceConfig &config) {
  m_config = config;
  m_initialized = true;
  return true;
}

/**
 * @brief Run stub inference (no-op).
 *
 * Sets all output values to zero. Useful for testing the data pipeline
 * without actual model inference.
 *
 * @param inputs  Input data (unused)
 * @param outputs Output data, will be zeroed
 * @param batch_size Number of samples in batch
 * @return true if initialized, false otherwise
 */
bool StubBackend::infer(const double *inputs, double *outputs, int batch_size) {
  (void)inputs; // Unused parameter

  if (!m_initialized) {
    return false;
  }

  // Zero all outputs
  const int output_size = batch_size * m_config.output_channels;
  std::memset(outputs, 0, output_size * sizeof(double));

  return true;
}

/**
 * @brief Finalize the stub backend.
 *
 * Resets initialization state. Backend can be re-initialized after this.
 */
void StubBackend::finalize() { m_initialized = false; }

} // namespace inference
} // namespace emulator
