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
 * Stores configuration for later use. Always succeeds for stub.
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
 * Leaves outputs unchanged. Useful for testing the data pipeline
 * without actual model inference.
 *
 * @param inputs  Input data (unused)
 * @param outputs Output data (unchanged)
 * @param batch_size Number of samples in batch (unused)
 * @return true if initialized, false otherwise
 */
bool StubBackend::infer(const double *inputs, double *outputs, int batch_size) {
  (void)inputs;     // Unused
  (void)outputs;    // Unchanged
  (void)batch_size; // Unused

  if (!m_initialized) {
    return false;
  }

  // No-op: leave outputs unchanged
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
