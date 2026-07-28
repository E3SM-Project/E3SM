/**
 * @file inference_backend.cpp
 * @brief Lifecycle and the flat-array convenience path.
 */

#include "inference_backend.hpp"

#include "inference_error.hpp"

#include <iostream>

namespace emulator {
namespace inference {

void InferenceBackend::initialize() {
  if (m_initialized) {
    return;
  }
  init_impl();
  m_initialized = true;
  if (m_config.verbose && m_context.is_root()) {
    std::cout << "[emulator::inference] " << name()
              << " ready: " << m_context.to_string() << "\n";
  }
}

bool InferenceBackend::infer(const TensorMap &inputs, TensorMap &outputs) {
  EMULATOR_INFER_REQUIRE(m_initialized,
                         name() << " backend used before initialize().");
  return infer_impl(inputs, outputs);
}

bool InferenceBackend::infer(const double *inputs, double *outputs,
                             int batch_size) {
  EMULATOR_INFER_REQUIRE(batch_size >= 0,
                         "Negative batch size " << batch_size << ".");
  EMULATOR_INFER_REQUIRE(
      m_config.input_channels > 0 && m_config.output_channels > 0,
      "The flat-array infer() needs input_channels and output_channels; got "
          << m_config.input_channels << " and " << m_config.output_channels
          << ". Use the tensor overload for models with several fields.");

  TensorMap in;
  in.wrap(m_config.inputs.empty() ? "input" : m_config.inputs.front(), inputs,
          {batch_size, m_config.input_channels});
  TensorMap out;
  out.wrap(m_config.outputs.empty() ? "output" : m_config.outputs.front(),
           outputs, {batch_size, m_config.output_channels});
  return infer(in, out);
}

void InferenceBackend::finalize() {
  if (!m_initialized) {
    return;
  }
  final_impl();
  m_initialized = false;
}

} // namespace inference
} // namespace emulator
