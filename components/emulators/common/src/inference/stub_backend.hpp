/**
 * @file stub_backend.hpp
 * @brief Stub inference backend for testing without ML dependencies.
 */

#ifndef E3SM_STUB_BACKEND_HPP
#define E3SM_STUB_BACKEND_HPP

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

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

  /// @copydoc InferenceBackend::infer
  bool infer(const double *inputs, double *outputs,
             int batch_size = 1) override;

  /// @copydoc InferenceBackend::finalize
  void finalize() override;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "Stub"; }
};

} // namespace inference
} // namespace emulator

#endif // E3SM_STUB_BACKEND_HPP
