/**
 * @file libtorch_backend.cpp
 * @brief LibTorch inference backend implementation.
 *
 * Provides native C++ neural network inference using LibTorch (PyTorch C++
 * API). This backend loads TorchScript models and executes inference without
 * Python.
 *
 * @note Models must be exported to TorchScript format (.pt) using
 * torch.jit.trace() or torch.jit.script() before use with this backend.
 *
 * @see LibTorchBackend
 */

#include "libtorch_backend.hpp"

// Standard library
#include <iostream>

// LibTorch headers
#include <torch/script.h>
#include <torch/torch.h>

namespace emulator {
namespace inference {

/**
 * @brief Private implementation details for LibTorchBackend.
 *
 * Uses PIMPL idiom to hide LibTorch types from the header, avoiding
 * the need to include torch headers in dependent code.
 */
struct LibTorchBackend::Impl {
  torch::jit::script::Module model;          ///< Loaded TorchScript model
  torch::Device device = torch::kCPU;        ///< Execution device (CPU or CUDA)
  torch::ScalarType dtype = torch::kFloat32; ///< Model precision
  bool model_loaded = false; ///< Whether model was successfully loaded
};

LibTorchBackend::LibTorchBackend() = default;

LibTorchBackend::~LibTorchBackend() {
  if (m_initialized) {
    finalize();
  }
}

/**
 * @brief Initialize the LibTorch backend.
 *
 * Loads the TorchScript model from the configured path and sets up
 * the execution device and precision.
 */
bool LibTorchBackend::initialize(const InferenceConfig &config) {
  m_config = config;
  m_impl = std::make_unique<Impl>();

  // Configure execution device
  if (m_config.device_id >= 0 && torch::cuda::is_available()) {
    m_impl->device = torch::Device(torch::kCUDA, m_config.device_id);
  } else {
    m_impl->device = torch::kCPU;
  }

  // Configure precision
  // Default to Float64 for E3SM compatibility (E3SM uses double precision)
  // FP16 or FP32 only available/recommended on CUDA for acceleration
  if (m_config.use_fp16 && m_impl->device.is_cuda()) {
    m_impl->dtype = torch::kFloat16;
  } else {
    // Use Float32 by default
    // NOTE: E3SM uses double precision by default
    // TODO: make this configurable if needed
    m_impl->dtype = torch::kFloat32;
  }

  // Load TorchScript model
  try {
    m_impl->model = torch::jit::load(m_config.model_path, m_impl->device);
    m_impl->model.eval(); // Set to evaluation mode (disables dropout, etc.)
  } catch (const c10::Error &e) {
    std::cerr << "[LibTorchBackend] Failed to load model: " << e.what()
              << std::endl;
    return false;
  }
  m_impl->model_loaded = true;
  m_initialized = true;
  return true;
}

/**
 * @brief Run inference on input data.
 *
 * Executes the loaded model on the provided input data. The backend
 * expects input in [batch_size, input_channels] format. For CNN models
 * requiring [N, C, H, W], the caller (EmulatorComp) must reshape first.
 *
 * @param inputs  Input data array of size [batch_size * input_channels]
 * @param outputs Output data array of size [batch_size * output_channels]
 * @param batch_size Number of samples in the batch
 * @return true if inference succeeded, false on error
 */
bool LibTorchBackend::infer(const double *inputs, double *outputs,
                            int batch_size) {
  if (!m_initialized || !m_impl || !m_impl->model_loaded) {
    std::cerr << "[LibTorchBackend::infer] ERROR: Backend not initialized!"
              << std::endl;
    return false;
  }

  const int C_in = m_config.input_channels;
  const int C_out = m_config.output_channels;

  try {
    torch::NoGradGuard no_grad;

    // Create input tensor by wrapping the pointer directly
    // EmulatorAtm has already arranged data in the correct shape:
    // - Spatial mode: [batch_size, C, H, W] flattened
    // - Pointwise mode: [batch_size, C]
    std::vector<int64_t> input_shape;
    if (m_config.spatial_mode) {
      input_shape = {batch_size, C_in, m_config.grid_height,
                     m_config.grid_width};
    } else {
      input_shape = {batch_size, C_in};
    }

    // We assume inputs are in double precision (Float64) and on host
    // TODO: make this configurable if needed
    torch::Tensor input_tensor = torch::from_blob(
        const_cast<double *>(inputs), input_shape,
        torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU));

    // Convert to model dtype (Float32) and target device
    // TODO: make EmulatorComp sends us the right dtype
    input_tensor = input_tensor.to(m_impl->device, m_impl->dtype);

    // Execute model forward pass
    std::vector<torch::jit::IValue> model_inputs;
    model_inputs.push_back(input_tensor);

    torch::Tensor output_tensor =
        m_impl->model.forward(model_inputs).toTensor();

    // Convert output back to Float64 on CPU
    // TODO: make this configurable if needed
    output_tensor = output_tensor.to(torch::kCPU, torch::kFloat64);

    // Ensure contiguous and copy
    if (!output_tensor.is_contiguous()) {
      output_tensor = output_tensor.contiguous();
    }

    // Output shape matches what EmulatorComp expects:
    // - Spatial: [batch_size, C, H, W]
    // - Pointwise: [batch_size, C]
    const size_t output_size = output_tensor.numel();
    std::memcpy(outputs, output_tensor.data_ptr<double>(),
                output_size * sizeof(double));

  } catch (const c10::Error &e) {
    std::cerr << "[LibTorchBackend::infer] FATAL c10::Error: " << e.what()
              << std::endl;
    return false;
  } catch (const std::exception &e) {
    std::cerr << "[LibTorchBackend::infer] FATAL exception: " << e.what()
              << std::endl;
    return false;
  }

  return true;
}

void LibTorchBackend::finalize() {
  m_impl.reset();
  m_initialized = false;
}

size_t LibTorchBackend::get_memory_usage_bytes() const {
  return m_model_memory_bytes;
}

} // namespace inference
} // namespace emulator
