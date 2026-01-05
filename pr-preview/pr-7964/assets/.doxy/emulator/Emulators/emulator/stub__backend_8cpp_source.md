

# File stub\_backend.cpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**inference**](dir_93d8073dc3188aa38589a4c97b36101e.md) **>** [**stub\_backend.cpp**](stub__backend_8cpp.md)

[Go to the documentation of this file](stub__backend_8cpp.md)


```C++


#include "stub_backend.hpp"
#include <cstring>

namespace emulator {
namespace inference {

bool StubBackend::initialize(const InferenceConfig &config) {
  m_config = config;
  m_initialized = true;
  return true;
}

bool StubBackend::infer(const double *inputs, double *outputs, int batch_size) {
  (void)inputs; // Unused parameter

  if (!m_initialized) {
    return false;
  }

  // Zero all outputs
  const int output_size = batch_size * m_config.output_channels;
  std::memset(outputs, 0, output_size * sizeof(double));

  return true;
}

void StubBackend::finalize() { m_initialized = false; }

} // namespace inference
} // namespace emulator
```


