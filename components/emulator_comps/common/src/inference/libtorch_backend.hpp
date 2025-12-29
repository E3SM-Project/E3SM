/**
 * @file libtorch_backend.hpp
 * @brief LibTorch inference backend for native C++ PyTorch inference.
 */

#ifndef LIBTORCH_BACKEND_HPP
#define LIBTORCH_BACKEND_HPP

#include "inference_backend.hpp"
#include <memory>
#include <string>

namespace emulator {
namespace inference {

/**
 * @brief LibTorch backend for native C++ PyTorch inference.
 *
 * Uses LibTorch (PyTorch C++ API) to load and run TorchScript models.
 * This provides production-grade inference without Python dependencies.
 *
 * ## Advantages
 * 
 * - Native C++ performance with no interpreter overhead
 * - No Python dependency at runtime
 * - Full GPU support via CUDA
 * - Thread-safe (no GIL limitations)
 *
 * ## Limitations
 * 
 * - Requires model export to TorchScript format
 * - LibTorch library must be available at build time
 * - Model must be compatible with torch.jit.trace() or torch.jit.script()
 *
 * ## Data Format
 * Input and output tensors are expected in [batch_size, channels] format.
 * The backend does NOT perform spatial reshaping - callers are responsible
 * for providing data in the format expected by their model.
 *
 * ## Configuration
 * 
 * - `model_path`: Path to TorchScript model (.pt file)
 * - `device_id`: GPU device ID (-1 for CPU, 0+ for CUDA device)
 * - `use_fp16`: Use half precision (requires CUDA, may improve performance)
 * - `input_channels`: Number of input features per sample
 * - `output_channels`: Number of output features per sample
 *
 * @see InferenceBackend for the base interface
 * @see create_backend() for factory function
 */
class LibTorchBackend : public InferenceBackend {
public:
  LibTorchBackend();
  ~LibTorchBackend() override;

  /// @copydoc InferenceBackend::initialize
  bool initialize(const InferenceConfig &config) override;

  /// @copydoc InferenceBackend::infer
  bool infer(const double *inputs, double *outputs, int batch_size) override;

  /// @copydoc InferenceBackend::finalize
  void finalize() override;

  /// @copydoc InferenceBackend::name
  std::string name() const override { return "LibTorch"; }

  /// @copydoc InferenceBackend::is_initialized
  bool is_initialized() const override { return m_initialized; }

  /// @copydoc InferenceBackend::type
  BackendType type() const override { return BackendType::LIBTORCH; }

  /**
   * @brief Get approximate memory usage in bytes.
   * @return Estimated memory used by model and buffers
   */
  size_t get_memory_usage_bytes() const;

private:
  bool m_initialized = false;      ///< Initialization state
  InferenceConfig m_config;        ///< Stored configuration
  size_t m_model_memory_bytes = 0; ///< Cached memory usage estimate

  /// @brief Private implementation (PIMPL idiom)
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

} // namespace inference
} // namespace emulator

#endif // LIBTORCH_BACKEND_HPP
