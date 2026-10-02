/**
 * @file inference_backend.cpp
 * @brief Flat-array convenience path shared by all backends.
 */

#include "inference_backend.hpp"

#include "inference_error.hpp"

namespace emulator {
namespace inference {

bool InferenceBackend::infer(const double *inputs, double *outputs,
                             int batch_size) {
  EMULATOR_INFER_REQUIRE(
      m_config.input_channels > 0 && m_config.output_channels > 0,
      "The flat-array infer() needs input_channels and output_channels, got "
          << m_config.input_channels << " and " << m_config.output_channels
          << ".");

  TensorMap in;
  in.wrap("input", inputs, {batch_size, m_config.input_channels});
  TensorMap out;
  out.wrap("output", outputs, {batch_size, m_config.output_channels});
  return infer(in, out);
}

} // namespace inference
} // namespace emulator
