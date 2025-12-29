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
 * Provides a no-op implementation that returns zeros for all outputs.
 * Useful for:
 * 
 * - CI/CD testing without ML dependencies
 * - Debugging coupling and I/O logic
 * - Development when GPU or PyTorch is unavailable
 * - Validating data pipeline correctness
 *
 * ## Behavior
 * 
 * - initialize() always succeeds
 * - infer() zeros all output values
 * - No external dependencies required
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
  bool infer(const double *inputs, double *outputs, int batch_size) override;

  /// @copydoc InferenceBackend::finalize
  void finalize() override;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "Stub"; }

  /// @copydoc InferenceBackend::is_initialized
  bool is_initialized() const override { return m_initialized; }

  /// @copydoc InferenceBackend::type
  BackendType type() const override { return BackendType::STUB; }

private:
  bool m_initialized = false; ///< Initialization state
  InferenceConfig m_config;   ///< Stored configuration
};

} // namespace inference
} // namespace emulator

#endif // STUB_BACKEND_HPP
