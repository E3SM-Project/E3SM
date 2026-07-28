/**
 * @file stub_inference_backend.cpp
 * @brief Stub inference backend implementation.
 */

#include "stub_inference_backend.hpp"

#include <iostream>

namespace emulator {
namespace inference {

/**
 * @brief Run stub inference: leave the outputs exactly as they were.
 *
 * Not zero-filling is deliberate.  A component that forgets to write a
 * field then shows the previous step's values rather than a plausible-looking
 * zero, which is much easier to notice.
 */
bool StubBackend::infer_impl(const TensorMap &inputs, TensorMap &outputs) {
  ++m_calls;
  if (m_config.verbose && m_context.is_root()) {
    std::cout << "[emulator::inference] stub step " << m_calls << ": "
              << inputs.size() << " input(s), " << outputs.size()
              << " output(s), unchanged\n";
  }
  return true;
}

} // namespace inference
} // namespace emulator
