/**
 * @file stub_inference_backend.hpp
 * @brief Stub inference backend for testing without ML dependencies.
 */

#ifndef E3SM_EMULATOR_STUB_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_STUB_INFERENCE_BACKEND_HPP

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Backend that runs no model, for CI and component bring-up.
 *
 * Leaves outputs unchanged, so it exercises every code path around a model
 * -- configuration, the context, tensor wrapping, the component's pack and
 * unpack -- with nothing installed.
 *
 * @see InferenceBackend for the base interface
 */
class StubBackend : public InferenceBackend {
public:
  StubBackend(const InferenceConfig &config, const InferenceContext &context)
      : InferenceBackend(config, context) {}
  ~StubBackend() override = default;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "Stub"; }

protected:
  void init_impl() override {}
  bool infer_impl(const TensorMap &inputs, TensorMap &outputs) override;
  void final_impl() override {}

private:
  int m_calls = 0;
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_STUB_INFERENCE_BACKEND_HPP
