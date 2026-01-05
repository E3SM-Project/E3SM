

# File stub\_backend.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**inference**](dir_93d8073dc3188aa38589a4c97b36101e.md) **>** [**stub\_backend.hpp**](stub__backend_8hpp.md)

[Go to the documentation of this file](stub__backend_8hpp.md)


```C++


#ifndef STUB_BACKEND_HPP
#define STUB_BACKEND_HPP

#include "inference_backend.hpp"

namespace emulator {
namespace inference {

class StubBackend : public InferenceBackend {
public:
  StubBackend() = default;
  ~StubBackend() override = default;

  bool initialize(const InferenceConfig &config) override;

  bool infer(const double *inputs, double *outputs, int batch_size) override;

  void finalize() override;

  std::string name() const override { return "Stub"; }

  bool is_initialized() const override { return m_initialized; }

  BackendType type() const override { return BackendType::STUB; }

private:
  bool m_initialized = false; 
  InferenceConfig m_config;   
};

} // namespace inference
} // namespace emulator

#endif // STUB_BACKEND_HPP
```


