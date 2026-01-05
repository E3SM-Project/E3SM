

# File inference\_factory.cpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**inference**](dir_93d8073dc3188aa38589a4c97b36101e.md) **>** [**inference\_factory.cpp**](inference__factory_8cpp.md)

[Go to the documentation of this file](inference__factory_8cpp.md)


```C++


#include "inference_backend.hpp"
#include "stub_backend.hpp"

#ifdef EMULATOR_HAS_LIBTORCH
#include "libtorch_backend.hpp"
#endif

namespace emulator {
namespace inference {

std::unique_ptr<InferenceBackend> create_backend(BackendType type) {
  switch (type) {
  case BackendType::STUB:
    return std::make_unique<StubBackend>();

#ifdef EMULATOR_HAS_LIBTORCH
  case BackendType::LIBTORCH:
    return std::make_unique<LibTorchBackend>();
#endif

  default:
    // Fall back to stub for unknown or unavailable backends
    return std::make_unique<StubBackend>();
  }
}

std::unique_ptr<InferenceBackend>
create_backend(const InferenceConfig &config) {
  auto backend = create_backend(config.backend);

  if (backend && !backend->initialize(config)) {
    return nullptr;
  }

  return backend;
}

} // namespace inference
} // namespace emulator
```


