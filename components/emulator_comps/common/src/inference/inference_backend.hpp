/**
 * @file inference_backend.hpp
 * @brief Abstract interface for neural network inference backends.
 *
 * Defines the common interface for all inference backends, allowing
 * pluggable AI/ML frameworks for model execution.
 */

#ifndef INFERENCE_BACKEND_HPP
#define INFERENCE_BACKEND_HPP

#include <memory>
#include <string>
#include <vector>

// MPI type handling:
// Try to include mpi.h if available to get real definitions.
// If mpi.h is not available, provide stub definitions for non-MPI builds.
#if defined(__has_include) && __has_include(<mpi.h>)
#include <mpi.h>
#else
// No mpi.h available - provide stub definitions for non-MPI builds
typedef int MPI_Comm;
#define MPI_COMM_WORLD 0
#define MPI_COMM_NULL 0
#endif

namespace emulator {
namespace inference {

/**
 * @brief Enumeration of available inference backend types.
 *
 * ## Current Backends
 *
 * - **STUB**: No-op backend for testing. Always returns zeros. Useful for
 *   development and CI testing without ML dependencies.
 *
 * - **LIBTORCH**: Uses LibTorch (PyTorch's C++ API) for inference. Requires
 *   a traced TorchScript model (.pt file). Best for production use with
 *   PyTorch models. Supports CPU and GPU execution.
 *
 * ## Planned Backends
 *
 * - **PYTORCH**: Will embed a Python interpreter directly to run native
 *   PyTorch models without tracing/scripting. Useful for rapid prototyping
 *   and models that don't trace well. Requires Python + PyTorch installation.
 *
 * - **ONNX**: Will use ONNX Runtime for cross-framework inference. Allows
 *   running models exported from PyTorch, TensorFlow, JAX, etc. in a
 *   unified, optimized runtime. Good for production with diverse model sources.
 *
 * - **LAPIS**: Will provide Kokkos-based inference for tight integration
 *   with E3SM's GPU/accelerator infrastructure. Designed for models that
 *   need to share memory with other Kokkos-enabled components (e.g., EAMxx)
 *   without data movement overhead.
 */
enum class BackendType {
  STUB,     ///< No-op backend for testing (no ML dependencies)
  LIBTORCH, ///< LibTorch C++ backend (native libtorch inference)
  // Future backends (uncomment when implemented):
  // PYTORCH, ///< Python interpreter backend (native PyTorch, no tracing)
  // ONNX,    ///< ONNX Runtime backend (cross-framework, optimized)
  // LAPIS,   ///< Kokkos interop backend (GPU memory sharing with E3SM)
};

/**
 * @brief Convert a BackendType to its string representation.
 * @param type The backend type
 * @return String name of the backend
 */
inline std::string backend_type_to_string(BackendType type) {
  switch (type) {
  case BackendType::STUB:
    return "STUB";
  case BackendType::LIBTORCH:
    return "LIBTORCH";
  default:
    return "UNKNOWN";
  }
}

/**
 * @brief Parse a backend type from a string.
 *
 * Supports case-insensitive matching and common aliases.
 *
 * @param str String to parse ("libtorch", "torch", "stub", etc.)
 * @return Corresponding BackendType (defaults to STUB if unrecognized)
 */
inline BackendType parse_backend_type(const std::string &str) {
  if (str == "LIBTORCH" || str == "libtorch" || str == "torch") {
    return BackendType::LIBTORCH;
  }
  return BackendType::STUB;
}

/**
 * @brief Configuration options for inference backends.
 *
 * Contains all parameters needed to initialize an inference backend,
 * including model path, device selection, and tensor dimensions.
 */
struct InferenceConfig {
  BackendType backend = BackendType::STUB; ///< Backend implementation to use

  std::string model_path; ///< Path to model file (TorchScript .pt for LibTorch)

  int device_id = -1;    ///< GPU device ID (-1 for CPU)
  bool use_fp16 = false; ///< Use half precision (requires CUDA)
  bool verbose = false;  ///< Enable verbose output (for debugging)

  int input_channels = 44;  ///< Number of input features per sample
  int output_channels = 50; ///< Number of output features per sample

  // Spatial mode settings (for CNN models)
  bool spatial_mode =
      false;           ///< If true, reshape to [N, C, H, W] for CNN models
  int grid_height = 0; ///< Height dimension (for spatial_mode)
  int grid_width = 0;  ///< Width dimension (for spatial_mode)

  // Validation and dry-run
  bool dry_run = false; ///< If true, validate config and exit without running
  std::vector<std::string>
      expected_input_vars; ///< Expected input variable names
  std::vector<std::string>
      expected_output_vars; ///< Expected output variable names
};

/**
 * @brief Result of configuration validation.
 *
 * Returned by InferenceBackend::validate() to report validation status
 * and any errors found.
 */
struct ValidationResult {
  bool valid = true;                 ///< True if all checks passed
  std::vector<std::string> errors;   ///< List of error messages
  std::vector<std::string> warnings; ///< List of warning messages

  /// Add an error and mark as invalid
  void add_error(const std::string &msg) {
    valid = false;
    errors.push_back(msg);
  }

  /// Add a warning (does not affect validity)
  void add_warning(const std::string &msg) { warnings.push_back(msg); }

  /// Check if there are any warnings
  bool has_warnings() const { return !warnings.empty(); }
};

/**
 * @brief Abstract interface for inference backends.
 *
 * Provides a unified API for running neural network inference
 * regardless of the underlying ML framework. All backends:
 *
 * - Accept input as a flat double array [batch_size * input_channels]
 * - Produce output as a flat double array [batch_size * output_channels]
 * - Handle device placement and precision conversion internally
 *
 * ## Thread Safety
 * Implementations should document their thread-safety characteristics.
 * In general, a single backend instance should only be used from one thread.
 *
 * @see LibTorchBackend, StubBackend for concrete implementations
 */
class InferenceBackend {
public:
  virtual ~InferenceBackend() = default;

  /**
   * @brief Initialize the backend.
   *
   * Loads the model, allocates resources, and prepares for inference.
   * Must be called before infer().
   *
   * @param config Configuration options
   * @return true if initialization succeeded, false on error
   */
  virtual bool initialize(const InferenceConfig &config) = 0;

  /**
   * @brief Run inference on input data.
   *
   * Executes the model on the provided input batch and writes results
   * to the output buffer.
   *
   * @param inputs  Input data array, size = batch_size * input_channels
   * @param outputs Output data array, size = batch_size * output_channels
   * @param batch_size Number of samples in the batch
   * @return true if inference succeeded, false on error
   *
   * @pre initialize() must have been called successfully
   * @pre outputs must be pre-allocated with sufficient size
   */
  virtual bool infer(const double *inputs, double *outputs, int batch_size) = 0;

  /**
   * @brief Release resources and finalize the backend.
   *
   * After calling this, the backend is no longer usable until
   * initialize() is called again.
   */
  virtual void finalize() = 0;

  /**
   * @brief Get the human-readable name of this backend.
   * @return Backend name (e.g., "LibTorch", "Stub")
   */
  virtual std::string name() const = 0;

  /**
   * @brief Check if the backend is ready for inference.
   * @return true if initialized and ready
   */
  virtual bool is_initialized() const = 0;

  /**
   * @brief Validate configuration before running.
   *
   * Checks that the model file exists, dimensions match, device is available,
   * etc. Call this after initialize() to detect configuration errors early.
   *
   * @return ValidationResult with errors/warnings if any
   */
  virtual ValidationResult validate() const {
    // Default implementation: always valid
    return ValidationResult{};
  }

  /**
   * @brief Get the backend type enumeration.
   * @return BackendType value
   */
  virtual BackendType type() const = 0;
};

/**
 * @brief Factory function to create an inference backend by type.
 *
 * Creates an uninitialized backend instance. Call initialize() on
 * the returned object before use.
 *
 * @param type Backend type to create
 * @return Unique pointer to new backend instance, or nullptr on error
 */
std::unique_ptr<InferenceBackend> create_backend(BackendType type);

/**
 * @brief Factory function to create and initialize a backend.
 *
 * Convenience function that creates a backend and calls initialize().
 *
 * @param config Configuration including backend type and init parameters
 * @return Initialized backend, or nullptr if creation/init failed
 */
std::unique_ptr<InferenceBackend> create_backend(const InferenceConfig &config);

} // namespace inference
} // namespace emulator

#endif // INFERENCE_BACKEND_HPP
