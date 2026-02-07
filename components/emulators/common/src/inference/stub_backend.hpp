/**
 * @file stub_backend.hpp
 * @brief Stub inference backend for testing without ML dependencies.
 */

#ifndef STUB_BACKEND_HPP
#define STUB_BACKEND_HPP

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Stub backend for testing without actual inference.
 *
 * @see InferenceBackend for the base interface
 */
class StubBackend : public InferenceBackend {
public:
  StubBackend() = default;
  ~StubBackend() override = default;

  /// @copydoc InferenceBackend::initialize
  bool initialize(const InferenceConfig &config) override;

  /// @copydoc InferenceBackend::infer
  bool infer(const double *inputs, double *outputs,
             int batch_size = 1) override;

  /// @copydoc InferenceBackend::finalize
  void finalize() override;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "Stub"; }

  /// @copydoc InferenceBackend::is_initialized
  bool is_initialized() const override { return m_initialized; }

private:
  bool m_initialized = false; ///< Initialization state
  InferenceConfig m_config;   ///< Stored configuration
};

} // namespace inference
} // namespace emulator

#endif // STUB_BACKEND_HPP
