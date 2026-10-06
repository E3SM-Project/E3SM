/**
 * @file libtorch_inference_backend.cpp
 * @brief LibTorch inference backend implementation.
 */

#include "libtorch_inference_backend.hpp"

#include "fpe_guard.hpp"
#include "inference_error.hpp"

#include <ATen/Parallel.h>
#include <torch/csrc/jit/python/update_graph_executor_opt.h>
#include <torch/cuda.h>
#include <torch/script.h>

#include <iostream>
#include <vector>

namespace emulator {
namespace inference {

namespace {

torch::Device parse_device(const std::string &spec) {
  torch::Device device(torch::kCPU);
  bool known = true;
  try {
    device = torch::Device(spec);
  } catch (const std::exception &) {
    known = false;
  }
  EMULATOR_INFER_REQUIRE(known && (device.is_cpu() || device.is_cuda()),
                         "Unknown device '" << spec
                                            << "'. Use cpu, cuda or cuda:N.");
  if (device.is_cuda()) {
    // Falling back to the CPU would hide a large slowdown.
    EMULATOR_INFER_REQUIRE(torch::cuda::is_available(),
                           "Device '" << spec
                                      << "' was requested but no CUDA device "
                                         "is available.");
    EMULATOR_INFER_REQUIRE(
        !device.has_index() ||
            device.index() < static_cast<int>(torch::cuda::device_count()),
        "Device '" << spec << "' was requested but only "
                   << torch::cuda::device_count() << " CUDA device(s) exist.");
  }
  return device;
}

torch::ScalarType parse_dtype(const std::string &spec) {
  if (spec == "float32") {
    return torch::kFloat32;
  }
  EMULATOR_INFER_REQUIRE(spec == "float64",
                         "Unknown dtype '" << spec
                                           << "'. Use float32 or float64.");
  return torch::kFloat64;
}

/// A torch tensor of the module's dtype, on the module's device.
at::Tensor to_torch(const Tensor &tensor, const torch::Device &device,
                    torch::ScalarType dtype) {
  const auto opts = torch::TensorOptions().dtype(torch::kFloat64);
  if (tensor.size() == 0) {
    return torch::empty(tensor.dims(), opts).to(device, dtype);
  }
  // A view, not a copy; the const_cast is safe because it is only read.
  // to() copies only if the device or dtype differ.
  return torch::from_blob(const_cast<double *>(tensor.cdata()), tensor.dims(),
                          opts)
      .to(device, dtype);
}

/// Copy a module output into the caller's tensor.
void from_torch(const at::Tensor &result, Tensor &tensor) {
  EMULATOR_INFER_REQUIRE(result.sizes().vec() == tensor.dims(),
                         "The module returned shape "
                             << result.sizes() << " for output "
                             << tensor.to_string() << ".");
  if (tensor.size() == 0) {
    return;
  }
  // copy_ handles the device transfer and the conversion back to double.
  torch::from_blob(tensor.data(), tensor.dims(),
                   torch::TensorOptions().dtype(torch::kFloat64))
      .copy_(result);
}

/// Flatten what forward() returned: a tensor, or a tuple or list of tensors.
std::vector<at::Tensor> collect_outputs(const torch::jit::IValue &result) {
  std::vector<at::Tensor> tensors;
  if (result.isTensor()) {
    tensors.push_back(result.toTensor());
  } else if (result.isTensorList()) {
    tensors = result.toTensorVector();
  } else if (result.isTuple()) {
    for (const auto &element : result.toTupleRef().elements()) {
      EMULATOR_INFER_REQUIRE(element.isTensor(),
                             "The module returned a tuple containing a "
                                 << element.tagKind() << ", not a tensor.");
      tensors.push_back(element.toTensor());
    }
  } else {
    EMULATOR_INFER_REQUIRE(false, "The module returned a "
                                      << result.tagKind()
                                      << "; expected a tensor, or a tuple or "
                                         "list of tensors.");
  }
  return tensors;
}

} // namespace

struct LibTorchBackend::Impl {
  torch::jit::Module module;
  torch::Device device{torch::kCPU};
  torch::ScalarType dtype = torch::kFloat32;
};

LibTorchBackend::LibTorchBackend(const InferenceConfig &config)
    : InferenceBackend(config), m_impl(new Impl()) {
  m_impl->device = parse_device(m_config.get("device", "cpu"));
  m_impl->dtype = parse_dtype(m_config.get("dtype", "float32"));

  const int num_threads = m_config.get_int("num_threads", 0);
  EMULATOR_INFER_REQUIRE(num_threads >= 0,
                         "num_threads must be non-negative, got "
                             << num_threads << ".");
  if (num_threads > 0) {
    at::set_num_threads(num_threads);
  }

  torch::jit::setGraphExecutorOptimize(
      m_config.get_bool("jit_optimize", false));

  EMULATOR_INFER_REQUIRE(!m_config.model_path.empty(),
                         "The LibTorch backend needs model_path.");
  try {
    FpeGuard no_fpe;
    m_impl->module = torch::jit::load(m_config.model_path, m_impl->device);
  } catch (const std::exception &e) {
    throw InferenceError("Could not load the TorchScript module '" +
                         m_config.model_path + "':\n" + e.what());
  }
  m_impl->module.eval();

  if (m_config.verbose) {
    std::cout << "[emulator::inference] LibTorch loaded "
              << m_config.model_path << " on " << m_impl->device << " as "
              << m_config.get("dtype", "float32") << "\n";
  }
}

LibTorchBackend::~LibTorchBackend() = default;

bool LibTorchBackend::infer(const TensorMap &inputs, TensorMap &outputs) {
  EMULATOR_INFER_REQUIRE(m_impl, "The LibTorch backend was finalized.");
  EMULATOR_INFER_REQUIRE(inputs.size() > 0,
                         "The LibTorch backend was given no input tensors.");

  torch::NoGradGuard no_grad;
  FpeGuard no_fpe;

  std::vector<torch::jit::IValue> args;
  args.reserve(inputs.size());
  for (const auto &tensor : inputs) {
    args.emplace_back(to_torch(tensor, m_impl->device, m_impl->dtype));
  }

  torch::jit::IValue result;
  try {
    result = m_impl->module.forward(args);
  } catch (const std::exception &e) {
    std::string shapes;
    for (const auto &tensor : inputs) {
      shapes += " " + tensor.to_string();
    }
    throw InferenceError("The TorchScript module '" + m_config.model_path +
                         "' failed in forward() with inputs" + shapes + ":\n" +
                         e.what());
  }

  const std::vector<at::Tensor> results = collect_outputs(result);
  EMULATOR_INFER_REQUIRE(results.size() == outputs.size(),
                         "The module returned " << results.size()
                                                << " tensor(s) but "
                                                << outputs.size()
                                                << " output(s) were given.");
  std::size_t i = 0;
  for (auto &tensor : outputs) {
    from_torch(results[i++], tensor);
  }
  return true;
}

void LibTorchBackend::finalize() { m_impl.reset(); }

} // namespace inference
} // namespace emulator
