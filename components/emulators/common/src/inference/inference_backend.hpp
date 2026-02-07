/**
 * @file inference_backend.hpp
 * @brief Abstract interface for neural network inference backends.
 */

#ifndef INFERENCE_BACKEND_HPP
#define INFERENCE_BACKEND_HPP

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
 */
class InferenceBackend {
public:
  virtual ~InferenceBackend() = default;

  /**
   * @brief Initialize the backend with configuration.
   * @param config Configuration options
   * @return true if initialization succeeded
   */
  virtual bool initialize(const InferenceConfig &config) = 0;

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

  /**
   * @brief Check if the backend is ready for inference.
   */
  virtual bool is_initialized() const = 0;
};

} // namespace inference
} // namespace emulator

#endif // INFERENCE_BACKEND_HPP
