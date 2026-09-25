/**
 * @file libtorch_inference_backend.hpp
 * @brief Inference backend that runs a TorchScript module through LibTorch.
 */

#ifndef E3SM_EMULATOR_LIBTORCH_INFERENCE_BACKEND_HPP
#define E3SM_EMULATOR_LIBTORCH_INFERENCE_BACKEND_HPP

#include <memory>

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

/**
 * @brief Runs a TorchScript module (saved by torch.jit.save) in-process.
 *
 * Input tensors are passed to forward() positionally, in TensorMap order,
 * with the shapes they declare. The module may return a tensor, or a tuple
 * or list of tensors; these are copied into the output tensors in order,
 * and their shapes must match exactly (nothing is reshaped to fit).
 *
 * Options read from `config.options`:
 * - `device`       cpu (default), cuda or cuda:N. An unavailable CUDA device
 *                  is an error, never a silent fallback to the CPU.
 * - `dtype`        float32 (default) or float64: the module's precision.
 *                  Tensors are converted on the way in and out.
 * - `num_threads`  Intra-op threads; 0 (default) keeps torch's default.
 * - `jit_optimize` Let TorchScript re-optimize the graph after the first
 *                  calls (default false). Faster, but results then change
 *                  between early and later calls, which breaks exact
 *                  restarts. Process-wide.
 *
 * @see InferenceBackend for the base interface
 */
class LibTorchBackend : public InferenceBackend {
public:
  /// @throws InferenceError if the options are invalid or the load fails
  explicit LibTorchBackend(const InferenceConfig &config);
  ~LibTorchBackend() override;

  using InferenceBackend::infer;

  /// @copydoc InferenceBackend::infer
  bool infer(const TensorMap &inputs, TensorMap &outputs) override;

  /// @copydoc InferenceBackend::finalize
  void finalize() override;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "LibTorch"; }

private:
  struct Impl; ///< Keeps torch headers out of this one
  std::unique_ptr<Impl> m_impl;
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_LIBTORCH_INFERENCE_BACKEND_HPP
