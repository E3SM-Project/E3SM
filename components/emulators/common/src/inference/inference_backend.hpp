/**
 * @file inference_backend.hpp
 * @brief Abstract interface for neural network inference backends.
 */

#ifndef E3SM_EMULATOR_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_INFERENCE_BACKEND_HPP

#include <string>

#include "inference_config.hpp"
#include "inference_context.hpp"
#include "tensor.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Abstract interface for inference backends.
 *
 * A backend evaluates a model given named input tensors and named output
 * tensors.  Whether that involves MPI collectives is a property of the
 * *model*, not of this interface: a column-local network needs none, while
 * a global model shards one sample across the component's ranks and its
 * `infer` is therefore collective -- every rank must call it the same
 * number of times, in the same order.  The context is what tells the model
 * which of those two worlds it is in.
 *
 * Lifecycle: construct, initialize() once, infer() per step, finalize().
 * Both ends are idempotent, because a component may finalize explicitly and
 * then be destroyed.
 */
class InferenceBackend {
public:
  InferenceBackend(const InferenceConfig &config,
                   const InferenceContext &context)
      : m_config(config), m_context(context) {}
  virtual ~InferenceBackend() = default;

  /// @brief Load the model.
  void initialize();

  /**
   * @brief Evaluate the model.
   * @param inputs  Named input tensors (usually const views of E3SM memory)
   * @param outputs Named output tensors (usually writable views of E3SM
   *                memory); written in place
   * @return true on success
   */
  bool infer(const TensorMap &inputs, TensorMap &outputs);

  /**
   * @brief Flat-array convenience overload, for code that thinks in arrays.
   *
   * Wraps both buffers as `[batch_size, channels]` tensors named by
   * `config.inputs[0]` / `config.outputs[0]` (defaulting to "input" and
   * "output") using `input_channels` / `output_channels`.
   */
  bool infer(const double *inputs, double *outputs, int batch_size = 1);

  /// @brief Release the model and any resources.
  void finalize();

  /// @brief Human-readable name of this backend.
  virtual std::string name() const = 0;

  bool is_initialized() const { return m_initialized; }
  const InferenceContext &context() const { return m_context; }

protected:
  /// @brief Load the model.  Called once, from initialize().
  virtual void init_impl() = 0;

  /// @brief Evaluate the model.  Called from infer(), after validation.
  virtual bool infer_impl(const TensorMap &inputs, TensorMap &outputs) = 0;

  /// @brief Release the model.  Called once, from finalize().
  virtual void final_impl() = 0;

  InferenceConfig m_config;   ///< Backend configuration
  InferenceContext m_context; ///< Ranks and decomposition from the coupler
  bool m_initialized = false;
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_INFERENCE_BACKEND_HPP
