# Inference Backends

The emulator framework uses a pluggable backend system for AI model inference.
This allows flexibility in deployment scenarios—from testing with no
dependencies to production with optimized C++ inference.

## Backend Interface

All backends implement the `InferenceBackend` interface:

```cpp
class InferenceBackend {
public:
    virtual bool initialize(const InferenceConfig& config) = 0;
    virtual bool infer(const double* inputs, double* outputs,
                       int batch_size) = 0;
    virtual void finalize() = 0;
    virtual std::string name() const = 0;
    virtual bool is_initialized() const = 0;
    virtual BackendType type() const = 0;
};
```

## Available Backends

| Backend | Purpose | Dependencies |
| ---------- | -------------------- | -------------------- |
| **STUB** | Testing without ML | None |
| **LibTorch** | Production inference | LibTorch C++ library |

---

## STUB Backend

**Purpose:** Testing and development without AI dependencies.

**Behavior:** Returns zeros for all outputs. Useful for validating the data
pipeline without actual model inference.

```cpp
InferenceConfig config;
config.backend = BackendType::STUB;
config.input_channels = 44;
config.output_channels = 50;

auto backend = create_backend(config);
// backend is already initialized by factory
```

**Use Cases:**

- CI/CD testing
- Debugging coupling logic
- Development without GPU or ML libraries

---

## LibTorch Backend

**Purpose:** Production-grade native C++ inference.

**How it works:**

1. Loads TorchScript model (`.pt` file) at initialization
2. Uses LibTorch C++ API for forward pass
3. No Python dependency at runtime

```cpp
InferenceConfig config;
config.backend = BackendType::LIBTORCH;
config.model_path = "/path/to/model.pt";
config.device_id = 0;     // GPU 0, or -1 for CPU
config.use_fp16 = false;  // Half precision (GPU only)
config.input_channels = 44;
config.output_channels = 50;

auto backend = create_backend(config);
```

### Data Format

The LibTorch backend operates on flat arrays with shape `[batch_size, channels]`:

- **Input:** `double[batch_size * input_channels]`
- **Output:** `double[batch_size * output_channels]`

    The backend does **not** perform spatial reshaping. For CNN models that
    expect `[N, C, H, W]` input, the caller (e.g., `EmulatorAtm`) must reshape
    the data appropriately before/after calling `infer()`.

### Spatial Mode (CNN Models)

For CNN-based models like ACE2, set `spatial_mode = true` in the emulator configuration:

```yaml
model_io:
  spatial_mode: true
  input_variables: [...]
  output_variables: [...]
```

When spatial mode is enabled, `EmulatorAtm` will:

1. Pack input fields from `[H*W, C]` to `[C, H, W]` (flattened as `[C*H*W]`)
2. Call backend with `batch_size=1` and `channels=C*H*W`
3. Unpack output from `[C, H, W]` back to `[H*W, C]`

### Preparing TorchScript Models

Export your PyTorch model to TorchScript format:

```python
import torch

model = YourModel()
model.load_state_dict(torch.load("weights.pt"))
model.eval()

# Option 1: Trace (for models with fixed control flow)
example_input = torch.randn(1, 44, 180, 360)
traced = torch.jit.trace(model, example_input)
traced.save("model_traced.pt")

# Option 2: Script (for models with dynamic control flow)
scripted = torch.jit.script(model)
scripted.save("model_scripted.pt")
```

### Build Requirements

Enable LibTorch support during CMake configuration:

```console
cmake -DEATM_INFERENCE_BACKEND=libtorch \
      -DTorch_DIR=/path/to/libtorch/share/cmake/Torch \
      ..
```

---

## Configuration Reference

| Option | Type | Default | Description |
| ---------------- | ------------ | ------- | -------------------------- |
| `backend` | `BackendType` | STUB | Backend implementation |
| `model_path` | string | "" | Path to TorchScript model |
| `device_id` | int | -1 | GPU device (-1 = CPU) |
| `use_fp16` | bool | false | Half precision (CUDA only) |
| `input_channels` | int | 44 | Input features per sample |
| `output_channels` | int | 50 | Output features per sample |
| `verbose` | bool | false | Enable debug output |

---

## Adding a New Backend

**Create header/source** in `common/src/inference/`:

```cpp
// my_backend.hpp
class MyBackend : public InferenceBackend {
public:
    bool initialize(const InferenceConfig& config) override;
    bool infer(const double* inputs, double* outputs, int batch_size) override;
    void finalize() override;
    
    std::string name() const override { return "MyBackend"; }
    bool is_initialized() const override { return m_initialized; }
    BackendType type() const override { return BackendType::MY_BACKEND; }
    
private:
    bool m_initialized = false;
    InferenceConfig m_config;
};
```

**Add to BackendType enum** in `inference_backend.hpp`:

```cpp
enum class BackendType {
    STUB,
    LIBTORCH,
    MY_BACKEND,  // New backend
};
```

**Register in factory** (`inference_factory.cpp`):

```cpp
case BackendType::MY_BACKEND:
    return std::make_unique<MyBackend>();
```

**Update CMakeLists.txt** to include new sources (conditionally if external
deps required)

**Add tests** in `tests/unit/inference/`
