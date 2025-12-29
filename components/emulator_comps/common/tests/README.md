# Inference Backend Testing Guide

This directory contains comprehensive tests for the E3SM emulator inference backends.

## Test Organization

```
tests/
├── unit/                          # Unit tests (single-rank)
│   ├── test_stub_backend.cpp     # Stub backend baseline
│   ├── test_inference_factory.cpp # Factory and type conversions
│   └── inference/                # Backend-specific tests
│       ├── test_libtorch_backend.cpp
│       └── test_pytorch_backend.cpp
└── integration/                   # Integration tests (multi-rank)
    └── test_distributed_backend.cpp # MPI distributed inference
```

## Quick Start (NERSC Perlmutter with ACE venv)

```bash
cd components/emulator_comps/common

# Run all tests (activates ACE venv automatically)
./run_inference_tests.sh

# Or run specific test
./run_inference_tests.sh test_pytorch_backend
```

The script will:
1. Activate ACE virtual environment (`/pscratch/sd/m/mahf708/ace/.venv`)
2. Generate dummy models
3. Build tests
4. Run all inference tests

## Running Tests Manually

### Prerequisites

**For PyTorch tests:**
```bash
# Activate ACE virtual environment
source /pscratch/sd/m/mahf708/ace/.venv/bin/activate

# Generate dummy models
cd tests/test_support/data
python3 create_dummy_model.py
cd ../../..
```

### On NERSC Perlmutter
```bash
# Interactive node
salloc -N 1 -C cpu --qos interactive -t 30:00

# Build in build directory
cd build
cmake ../components/emulator_comps/common \
  -DCMAKE_BUILD_TYPE=Debug \
  -DEMULATOR_HAS_MPI=ON

make
ctest --output-on-failure

# Run distributed test
srun -n 8 ./tests/test_distributed_backend
```

## Test Coverage

### ✅ Implemented Tests

**Stub Backend** (`test_stub_backend.cpp`)
- Initialization and finalization
- Basic inference (returns zeros)
- Configuration handling

**Inference Factory** (`test_inference_factory.cpp`)
- Backend creation for all types
- String-to-type conversion
- Type-to-string conversion
- Config-based creation

**Distributed Backend** (`test_distributed_backend.cpp`)
- MPI initialization
- Batch distribution across ranks
- Result gathering
- Multi-rank coordination

###🔧 TODO: Advanced Tests

**LibTorch Backend** (requires LibTorch + dummy TorchScript model)
- Model loading
- CPU inference
- GPU inference (if available)
- Batch processing

**ONNX Backend** (requires ONNX Runtime + dummy ONNX model)
- Model loading
- CPU/GPU execution providers
- Dynamic batch sizes

**PyTorch Backend** (requires Python/PyTorch environment)
- Python interoperability
- Module loading
- Distributed torch.distributed integration

**Performance Tests**
- Throughput benchmarks
- Latency measurements
- Scaling efficiency (weak/strong)

## Creating Test Models

### Dummy TorchScript Model
```python
import torch

# Simple linear model
model = torch.nn.Linear(10, 5)
model.eval()

# Trace and save
dummy_input = torch.randn(1, 10)
scripted = torch.jit.trace(model, dummy_input)
scripted.save("dummy_model.pt")
```

### Dummy ONNX Model
```python
import torch
import torch.onnx

model = torch.nn.Linear(10, 5)
dummy_input = torch.randn(1, 10)

torch.onnx.export(
    model, dummy_input, "dummy_model.onnx",
    opset_version=14,
    input_names=['input'], output_names=['output'],
    dynamic_axes={'input': {0: 'batch'}, 'output': {0: 'batch'}}
)
```

## Debugging Failed Tests

### Enable Verbose Output
```bash
# CMake verbose
cmake --build . --verbose

# CTest verbose
ctest --verbose --output-on-failure

# Direct test execution
./tests/test_stub_backend --verbose
```

### MPI Debugging
```bash
# GDB with MPI
mpirun -n 4 xterm -e gdb ./tests/test_distributed_backend

# Valgrind with MPI
mpirun -n 2 valgrind --leak-check=full ./tests/test_distributed_backend
```

## CI/CD Integration

Tests are designed to work in continuous integration:

```yaml
# Example GitHub Actions
- name: Run Tests
  run: |
    cd build
    ctest --output-on-failure
```

Set `CI=true` environment variable to skip interactive tests.

## Adding New Tests

1. Create test file in `unit/` or `integration/`
2. Add to `CMakeLists.txt`:
   ```cmake
   add_executable(test_my_feature unit/test_my_feature.cpp)
   target_link_libraries(test_my_feature PRIVATE emulator_common)
   add_test(NAME my_feature COMMAND test_my_feature)
   ```
3. Use assertions and clear output:
   ```cpp
   assert(result == expected);
   std::cout << "✓ Test passed" << std::endl;
   ```

## Test Data

Place test data files in `test_support/data/`:
- Small dummy models (< 1MB)
- Sample input/output arrays
- Configuration files

Access via:
```cpp
#include "../test_support/test_config.hpp"
std::string path = emulator::testing::get_test_data_file("dummy_model.pt");
```
