/**
 * @file stub_backend.hpp
 * @brief Stub inference backend for testing without ML dependencies.
 */

#ifndef E3SM_EMULATOR_STUB_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_STUB_INFERENCE_BACKEND_HPP

#include "inference_backend.hpp"

namespace e3sm::emulator::inference {

/**
 * @brief Stub backend for testing without actual inference.
 *
 * Leaves outputs unchanged. Useful for testing the data pipeline
 * without actual model inference.
 *
 * @see InferenceBackend for the base interface
 */
class StubBackend : public InferenceBackend {
public:
  explicit StubBackend(const InferenceConfig &config)
      : InferenceBackend(config) {}
  ~StubBackend() override = default;

  using InferenceBackend::infer;

  /// @copydoc InferenceBackend::infer
  bool infer(const TensorMap &inputs, TensorMap &outputs) override;

  /// @copydoc InferenceBackend::finalize
  void finalize() override;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "Stub"; }
};

} // namespace e3sm::emulator::inference

#endif // E3SM_EMULATOR_STUB_INFERENCE_BACKEND_HPP
