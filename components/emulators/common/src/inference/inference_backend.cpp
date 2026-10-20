/**
 * @file inference_backend.cpp
 * @brief Flat-array convenience path shared by all backends.
 */

#include "inference_backend.hpp"

#include "inference_error.hpp"

namespace e3sm::emulator::inference {

bool InferenceBackend::infer(const double *inputs, double *outputs, int batch_size) {
  EMULATOR_INFER_REQUIRE(
      config.input_channels > 0 && config.output_channels > 0,
      "The flat-array infer() needs input_channels and output_channels, got "
          << config.input_channels << " and " << config.output_channels
          << ".");

  TensorMap in;
  in.wrap("input", inputs, {batch_size, config.input_channels});
  TensorMap out;
  out.wrap("output", outputs, {batch_size, config.output_channels});
  return infer(in, out);
}

} // namespace e3sm::emulator::inference
