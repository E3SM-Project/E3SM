/**
 * @file stub_backend.cpp
 * @brief Stub inference backend implementation.
 *
 * Provides a no-op backend for testing the emulator framework without
 * requiring ML dependencies. Leaves outputs unchanged.
 */

#include "stub_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Run stub inference (no-op).
 *
 * Leaves outputs unchanged. Useful for testing the data pipeline
 * without actual model inference.
 *
 * @param inputs  Input data (unused)
 * @param outputs Output data (unchanged)
 * @param batch_size Number of samples in batch (unused)
 * @return Always returns true
 */
bool StubBackend::infer(const double *inputs, double *outputs, int batch_size) {
  (void)inputs;     // Unused
  (void)outputs;    // Unchanged
  (void)batch_size; // Unused

  // No-op: leave outputs unchanged
  return true;
}

/**
 * @brief Finalize the stub backend.
 *
 * No-op for stub backend; included for interface compliance.
 */
void StubBackend::finalize() {}

} // namespace inference
} // namespace emulator
