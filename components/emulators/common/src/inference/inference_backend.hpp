/**
 * @file inference_backend.hpp
 * @brief Abstract interface for neural network inference backends.
 */

#ifndef E3SM_EMULATOR_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_INFERENCE_BACKEND_HPP

#include <string>

namespace emulator {
namespace inference {

/**
 * @brief Enumeration of available inference backend types.
 */
enum class BackendType {
  STUB, ///< No-op backend for testing (no ML dependencies)
};

/**
 * @brief Minimal configuration for inference backends.
 */
struct InferenceConfig {
  int input_channels = 0;  ///< Number of input features per grid point
  int output_channels = 0; ///< Number of output features per grid point
  bool verbose = false;    ///< Enable verbose output (for debugging)
};

/**
 * @brief Abstract interface for inference backends.
 *
 * Provides a unified API for running neural network inference.
 * All backends:
 * - Accept input as a flat array [batch_size * input_channels]
 * - Produce output as a flat array [batch_size * output_channels]
 *
 * Backends are fully configured on construction (config is passed
 * to the constructor). No separate initialization step is needed.
 */
class InferenceBackend {
public:
  explicit InferenceBackend(const InferenceConfig &config) : m_config(config) {}
  virtual ~InferenceBackend() = default;

  /**
   * @brief Run inference on input data.
   * @param inputs  Input data array [batch_size * input_channels]
   * @param outputs Output data array [batch_size * output_channels]
   * @param batch_size Number of samples in the batch
   * @return true if inference succeeded
   */
  virtual bool infer(const double *inputs, double *outputs,
                     int batch_size = 1) = 0;

  /**
   * @brief Release resources and finalize the backend.
   */
  virtual void finalize() = 0;

  /**
   * @brief Get the human-readable name of this backend.
   */
  virtual std::string name() const = 0;

protected:
  InferenceConfig m_config; ///< Backend configuration
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_INFERENCE_BACKEND_HPP
