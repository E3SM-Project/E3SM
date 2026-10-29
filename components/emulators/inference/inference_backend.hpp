/**
 * @file inference_backend.hpp
 * @brief Abstract interface for neural network inference backends.
 */

#ifndef E3SM_EMULATOR_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_INFERENCE_BACKEND_HPP

#include <memory>
#include <string>

#include "inference_config.hpp"
#include "tensor.hpp"

namespace e3sm::emulator::inference {

/**
 * @brief Enumeration of available inference backend types.
 */
enum class BackendType {
  STUB,     ///< No-op backend for testing (no ML dependencies)
  LIBTORCH, ///< TorchScript module run through LibTorch
  PYTHON,   ///< Python model run through an embedded interpreter
};

/**
 * @brief Abstract interface for inference backends.
 *
 * Provides a unified API for running neural network inference on named
 * tensors. Nothing here is specific to the caller: a component and a
 * single parameterization use it the same way.
 *
 * Backends are fully configured on construction (config is passed
 * to the constructor). No separate initialization step is needed.
 */
class InferenceBackend {
public:
  explicit InferenceBackend(const InferenceConfig &config) : config_(config) {}
  virtual ~InferenceBackend() = default;

  /**
   * @brief Access the configuration for this backend.
   */
  const InferenceConfig &config() const { return config_; }

  /**
   * @brief Run inference.
   * @param inputs  Input tensors (usually const views of model memory)
   * @param outputs Output tensors, written in place
   * @return true if inference succeeded
   * @throws InferenceError on misuse, e.g. a shape mismatch
   */
  virtual bool infer(const TensorMap &inputs, TensorMap &outputs) = 0;

  /**
   * @brief Run inference on flat arrays.
   *
   * Wraps the arrays as [batch_size, input_channels] and
   * [batch_size, output_channels] tensors and calls the overload above.
   *
   * @param inputs  Input data array [batch_size * input_channels]
   * @param outputs Output data array [batch_size * output_channels]
   * @param batch_size Number of samples in the batch
   * @return true if inference succeeded
   */
  bool infer(const double *inputs, double *outputs, int batch_size = 1);

  /**
   * @brief Release resources and finalize the backend.
   */
  virtual void finalize() = 0;

  /**
   * @brief Get the human-readable name of this backend.
   */
  virtual std::string name() const = 0;

private:
  InferenceConfig config_; ///< Backend configuration
};

} // namespace emulator::inference

#endif // E3SM_EMULATOR_INFERENCE_BACKEND_HPP
