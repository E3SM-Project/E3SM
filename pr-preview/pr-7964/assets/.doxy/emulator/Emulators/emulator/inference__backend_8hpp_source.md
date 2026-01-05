

# File inference\_backend.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**inference**](dir_93d8073dc3188aa38589a4c97b36101e.md) **>** [**inference\_backend.hpp**](inference__backend_8hpp.md)

[Go to the documentation of this file](inference__backend_8hpp.md)


```C++


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

enum class BackendType {
  STUB,     
  LIBTORCH, 
  // Future backends (uncomment when implemented):
  // PYTORCH, ///< Python interpreter backend (native PyTorch, no tracing)
  // ONNX,    ///< ONNX Runtime backend (cross-framework, optimized)
  // LAPIS,   ///< Kokkos interop backend (GPU memory sharing with E3SM)
};

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

inline BackendType parse_backend_type(const std::string &str) {
  if (str == "LIBTORCH" || str == "libtorch" || str == "torch") {
    return BackendType::LIBTORCH;
  }
  return BackendType::STUB;
}

struct InferenceConfig {
  BackendType backend = BackendType::STUB; 

  std::string model_path; 

  int device_id = -1;    
  bool use_fp16 = false; 
  bool verbose = false;  

  int input_channels = 44;  
  int output_channels = 50; 

  // Spatial mode settings (for CNN models)
  bool spatial_mode =
      false;           
  int grid_height = 0; 
  int grid_width = 0;  

  // Validation and dry-run
  bool dry_run = false; 
  std::vector<std::string>
      expected_input_vars; 
  std::vector<std::string>
      expected_output_vars; 
};

struct ValidationResult {
  bool valid = true;                 
  std::vector<std::string> errors;   
  std::vector<std::string> warnings; 

  void add_error(const std::string &msg) {
    valid = false;
    errors.push_back(msg);
  }

  void add_warning(const std::string &msg) { warnings.push_back(msg); }

  bool has_warnings() const { return !warnings.empty(); }
};

class InferenceBackend {
public:
  virtual ~InferenceBackend() = default;

  virtual bool initialize(const InferenceConfig &config) = 0;

  virtual bool infer(const double *inputs, double *outputs, int batch_size) = 0;

  virtual void finalize() = 0;

  virtual std::string name() const = 0;

  virtual bool is_initialized() const = 0;

  virtual ValidationResult validate() const {
    // Default implementation: always valid
    return ValidationResult{};
  }

  virtual BackendType type() const = 0;
};

std::unique_ptr<InferenceBackend> create_backend(BackendType type);

std::unique_ptr<InferenceBackend> create_backend(const InferenceConfig &config);

} // namespace inference
} // namespace emulator

#endif // INFERENCE_BACKEND_HPP
```


