/**
 * @file inference_backend.cpp
 * @brief Flat-array convenience path shared by all backends.
 */

#include "inference_backend.hpp"
#include "inference_error.hpp"
#include <memory>

namespace e3sm::emulator::inference {

bool InferenceBackend::infer(const double *inputs, double *outputs, int batch_size) {
  EMULATOR_INFER_REQUIRE(
      config_.input_channels > 0 && config_.output_channels > 0,
      "The flat-array infer() needs input_channels and output_channels, got "
          << config_.input_channels << " and " << config_.output_channels
          << ".");

  TensorMap in;
  in.wrap("input", inputs, {batch_size, config_.input_channels});
  TensorMap out;
  out.wrap("output", outputs, {batch_size, config_.output_channels});
  return infer(in, out);
}

} // namespace e3sm::emulator::inference
