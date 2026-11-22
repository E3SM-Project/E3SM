/**
 * @file libtorch_inference_backend.cpp
 * @brief LibTorch inference backend implementation.
 */

#include "libtorch_inference_backend.hpp"

#include "fpe_guard.hpp"
#include "inference_error.hpp"

#include <ATen/Parallel.h>
#include <ATen/core/ATen_fwd.h>
#include <c10/core/ScalarType.h>
#include <torch/csrc/jit/python/update_graph_executor_opt.h>
#include <torch/cuda.h>
#include <torch/script.h>

#include <algorithm>
#include <iostream>
#include <ranges>
#include <type_traits>
#include <vector>

namespace e3sm::emulator::inference {

namespace {

torch::Device parse_device(const std::string& spec) {
  torch::Device device(torch::kCPU);
  bool known = true;
  try {
    device = torch::Device(spec);
  } catch (const std::exception&) {
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

torch::ScalarType parse_dtype(const std::string& spec) {
  if (spec == "float32") {
    return torch::kFloat32;
  }
  EMULATOR_INFER_REQUIRE(spec == "float64",
                         "Unknown dtype '" << spec
                                           << "'. Use float32 or float64.");
  return torch::kFloat64;
}

/// A torch tensor of the module's dtype, on the module's device.
at::Tensor to_torch(const Tensor& tensor, const torch::Device& device) {

  return tensor.view().visit([&](auto values) -> at::Tensor {
    using T = typename decltype(values)::value_type;

    constexpr auto source_dtype = torch::CppTypeToScalarType<T>::value;

    // Note: opts device is torch::kCPU by default
    const auto opts = torch::TensorOptions().dtype(source_dtype);
    if (values.size() == 0) {
      return torch::empty(tensor.dims(), opts).to(device, source_dtype);
    }
    return torch::from_blob(const_cast<T*>(values.data()), tensor.dims(), opts)
        .to(device, source_dtype);
  });
}

/// Copy a module output into the caller's tensor.
void from_torch(const at::Tensor& result, Tensor& tensor) {

  EMULATOR_INFER_REQUIRE(std::ranges::equal(result.sizes(), tensor.dims()),
                         "The module returned shape "
                             << result.sizes() << " for output "
                             << tensor.to_string() << ".");
  if (tensor.size() == 0) {
    return;
  }
  tensor.view().visit([&](auto values) {
    using T = typename decltype(values)::value_type;
    constexpr auto dtype = torch::CppTypeToScalarType<T>::value;

    // copy_ handles the device transfer and the conversion back to the CPP
    // type.
    torch::from_blob(values.data(), tensor.dims(),
                     torch::TensorOptions().dtype(dtype))
        .copy_(result);
  });
}

/// Flatten what forward() returned: a tensor, or a tuple or list of tensors.
std::vector<at::Tensor> collect_outputs(const torch::jit::IValue& result) {
  std::vector<at::Tensor> tensors;
  if (result.isTensor()) {
    tensors.push_back(result.toTensor());
  } else if (result.isTensorList()) {
    tensors = result.toTensorVector();
  } else if (result.isTuple()) {
    for (const auto& element : result.toTupleRef().elements()) {
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

LibTorchBackend::LibTorchBackend(const InferenceConfig& config)
    : InferenceBackend(config), impl_(new Impl()) {
  impl_->device = parse_device(config.get("device", "cpu"));
  impl_->dtype = parse_dtype(config.get("dtype", "float32"));

  const int num_threads = config.get_int("num_threads", 0);
  EMULATOR_INFER_REQUIRE(num_threads >= 0,
                         "num_threads must be non-negative, got " << num_threads
                                                                  << ".");
  if (num_threads > 0) {
    at::set_num_threads(num_threads);
  }

  torch::jit::setGraphExecutorOptimize(config.get_bool("jit_optimize", false));

  EMULATOR_INFER_REQUIRE(!config.model_path.empty(),
                         "The LibTorch backend needs model_path.");
  try {
    FpeGuard no_fpe;
    impl_->module = torch::jit::load(config.model_path, impl_->device);
  } catch (const std::exception& e) {
    throw InferenceError("Could not load the TorchScript module '" +
                         config.model_path + "':\n" + e.what());
  }
  impl_->module.eval();

  if (config.verbose) {
    std::cout << "[emulator::inference] LibTorch loaded " << config.model_path
              << " on " << impl_->device << " as "
              << config.get("dtype", "float32") << "\n";
  }
}

LibTorchBackend::~LibTorchBackend() = default;

bool LibTorchBackend::infer(const TensorMap& inputs, TensorMap& outputs) {
  EMULATOR_INFER_REQUIRE(impl_, "The LibTorch backend was finalized.");
  EMULATOR_INFER_REQUIRE(inputs.size() > 0,
                         "The LibTorch backend was given no input tensors.");

  torch::NoGradGuard no_grad;
  FpeGuard no_fpe;

  std::vector<torch::jit::IValue> args;
  args.reserve(inputs.size());
  for (const auto& tensor : inputs) {
    args.emplace_back(to_torch(tensor, impl_->device));
  }

  torch::jit::IValue result;
  try {
    result = impl_->module.forward(args);
  } catch (const std::exception& e) {
    std::string shapes;
    for (const auto& tensor : inputs) {
      shapes += " " + tensor.to_string();
    }
    throw InferenceError("The TorchScript module '" + config().model_path +
                         "' failed in forward() with inputs" + shapes + ":\n" +
                         e.what());
  }

  const std::vector<at::Tensor> results = collect_outputs(result);
  EMULATOR_INFER_REQUIRE(results.size() == outputs.size(),
                         "The module returned "
                             << results.size() << " tensor(s) but "
                             << outputs.size() << " output(s) were given.");
  std::size_t i = 0;
  for (auto& tensor : outputs) {
    from_torch(results[i++], tensor);
  }
  return true;
}

void LibTorchBackend::finalize() { impl_.reset(); }

} // namespace e3sm::emulator::inference
