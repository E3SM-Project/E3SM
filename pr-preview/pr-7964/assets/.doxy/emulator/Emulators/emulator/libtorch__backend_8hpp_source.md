

# File libtorch\_backend.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**inference**](dir_93d8073dc3188aa38589a4c97b36101e.md) **>** [**libtorch\_backend.hpp**](libtorch__backend_8hpp.md)

[Go to the documentation of this file](libtorch__backend_8hpp.md)


```C++


#ifndef LIBTORCH_BACKEND_HPP
#define LIBTORCH_BACKEND_HPP

#include "inference_backend.hpp"
#include <memory>
#include <string>

namespace emulator {
namespace inference {

class LibTorchBackend : public InferenceBackend {
public:
  LibTorchBackend();
  ~LibTorchBackend() override;

  bool initialize(const InferenceConfig &config) override;

  bool infer(const double *inputs, double *outputs, int batch_size) override;

  void finalize() override;

  std::string name() const override { return "LibTorch"; }

  bool is_initialized() const override { return m_initialized; }

  BackendType type() const override { return BackendType::LIBTORCH; }

  size_t get_memory_usage_bytes() const;

private:
  bool m_initialized = false;      
  InferenceConfig m_config;        
  size_t m_model_memory_bytes = 0; 

  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

} // namespace inference
} // namespace emulator

#endif // LIBTORCH_BACKEND_HPP
```


