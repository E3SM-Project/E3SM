# Testing Guide

This guide covers the standalone testing framework for emulator_comps.

## Overview

The testing framework provides three levels of tests that can be run independently
of CIME, allowing for fast iteration during development and efficient CI.

| Level | Type | Description | Typical Runtime |
| ----- | ---- | ----------- | --------------- |
| 1 | Unit | Individual functions/classes | < 1 second |
| 2 | Component | Single component integration | 1-10 seconds |
| 3 | System | Full EATM end-to-end | 10-60 seconds |

## Quick Start

### Using the Build Script (Recommended)

```console
cd components/emulator_comps

# Build and run all tests
./test rebuild

# Run tests only
./test test

# Unit tests only
./test test -t1

# Quiet mode (minimal output)
./test rebuild -j8 -q
```

### Manual Build

```console
cd components/emulator_comps
mkdir build && cd build

# Configure with standalone mode and tests
cmake .. \
  -DEMULATOR_STANDALONE_BUILD=ON \
  -DEMULATOR_BUILD_TESTS=ON \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90

# Build
make -j8

# Run all tests
ctest --output-on-failure
```

### Running Specific Test Levels

```console
# Unit tests only (fastest, ~30 seconds)
ctest -L unit --output-on-failure

# Component tests
ctest -L component --output-on-failure

# System tests
ctest -L system --output-on-failure

# Tests for specific component
ctest -L eatm --output-on-failure

# MPI tests only (requires compute node)
ctest -L mpi --output-on-failure
```

> **Note:** MPI tests require `srun` which only works on compute nodes.
> Use `./test rebuild --mpi` in an salloc session.

### Test Level Configuration

Control which tests are built at configure time:

```console
# Only build unit tests (level 1)
cmake .. -DEMULATOR_TEST_LEVEL=1

# Build unit + component tests (levels 1-2)
cmake .. -DEMULATOR_TEST_LEVEL=2

# Build all tests (levels 1-3, default)
cmake .. -DEMULATOR_TEST_LEVEL=3
```

## Test Directory Structure

Tests are organized under each component:

```console
common/tests/
├── test_support/          # Shared test utilities
│   ├── emulator_test_session.cpp  # MPI-aware Catch2 main
│   ├── test_config.hpp    # Test paths configuration
│   └── test_data.hpp      # Synthetic data generators
├── unit/
│   ├── inference/         # Inference backend tests
│   │   ├── test_stub_backend.cpp
│   │   └── test_inference_factory.cpp
│   ├── test_emulator_config.cpp
│   └── test_coupling_fields.cpp
└── CMakeLists.txt

eatm/tests/
├── unit/                  # EATM-specific unit tests
│   ├── test_atm_field_manager.cpp
│   └── test_atm_coupling.cpp
├── component/             # Component integration tests
│   └── test_emulator_atm.cpp
├── system/                # Full system tests
│   ├── test_eatm_standalone.cpp
│   └── test_multi_timestep.cpp
└── CMakeLists.txt
```

## Writing Tests

### Using the Test Macros

The `EmulatorTestUtils.cmake` module provides convenient macros:

```cmake
# Unit test (level 1) - Fast, isolated tests
EmulatorUnitTest(test_my_function
  SOURCES test_my_function.cpp
  LIBS emulator_common
  LABELS myfeature
)

# Component test (level 2) - Single component
EmulatorComponentTest(test_my_component
  SOURCES test_my_component.cpp
  LIBS eatm
  LABELS eatm
)

# System test (level 3) with MPI
EmulatorSystemTest(test_full_run
  SOURCES test_full_run.cpp
  LIBS eatm
  LABELS eatm system
  MPI_RANKS 1 2 4
)
```

### Test Template

```cpp
#include <catch2/catch.hpp>
#include "test_data.hpp"

using namespace emulator::testing;

TEST_CASE("MyFeature basic functionality", "[unit][myfeature]") {
  
  SECTION("can do X") {
    // Setup
    auto grid = create_test_grid(100);
    
    // Exercise
    auto result = my_function(grid);
    
    // Verify
    REQUIRE(result.size() == 100);
  }
  
  SECTION("handles edge case Y") {
    REQUIRE_THROWS(my_function(nullptr));
  }
}
```

### MPI Tests

```cpp
#include <catch2/catch.hpp>
#include <mpi.h>
#include "test_data.hpp"

TEST_CASE("Parallel operation", "[component][mpi]") {
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  REQUIRE(nprocs >= 2);
  
  // Create partitioned grid
  auto grid = create_partitioned_grid(100, rank, nprocs);
  
  // Run parallel operation
  // ...
  
  // Collect and verify results
  int total_cols;
  MPI_Allreduce(&grid.ncol, &total_cols, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  REQUIRE(total_cols == 100);
}
```

## Test Data

### Synthetic Data

Use the `test_data.hpp` utilities for creating test data:

```cpp
#include "test_data.hpp"

using namespace emulator::testing;

// Create uniform grid
auto grid = create_test_grid(100, 72);

// Create partitioned grid for MPI
auto grid = create_partitioned_grid(1000, rank, nprocs);

// Generate random field (0 to 1)
auto field = create_random_field(100);

// Generate temperature-like field (250-320 K)
auto temp = create_temperature_field(100);

// Generate constant field
auto field = create_constant_field(100, 300.0);

// Compare arrays with tolerance
bool equal = arrays_equal(a.data(), b.data(), n, 1e-10, 1e-14);
```

### Real Data

For tests requiring real NetCDF data, place files in the component's `tests/data/`
directory. Access via:

```cpp
#include "test_config.hpp"

std::string path = get_test_data_file("test_grid.nc");
```

## CI Integration

### Test Labels

| Label | Description | Typical Count |
| ----- | ----------- | ------------- |
| `unit` | Unit tests | 10+ |
| `component` | Component tests | 5+ |
| `system` | System tests | 2-5 |
| `mpi` | Tests requiring multiple ranks | 5+ |
| `eatm` | EATM-related tests | 10+ |
| `inference` | Inference backend tests | 3+ |

### Example CI Workflows

```console
# Quick PR check (~30 seconds)
ctest -L unit --output-on-failure

# Standard CI (~2 minutes)
ctest -L "unit|component" --output-on-failure

# Full test suite (~5 minutes)
ctest --output-on-failure

# Nightly with all MPI configurations
for np in 1 2 4 8; do
  mpirun -np $np ctest -L mpi
done
```

### GitHub Actions Example

```yaml
test:
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - name: Build
      run: |
        mkdir build && cd build
        cmake .. -DEMULATOR_STANDALONE_BUILD=ON -DEMULATOR_BUILD_TESTS=ON
        make -j$(nproc)
    - name: Test
      run: |
        cd build
        ctest --output-on-failure -L unit
```

## Troubleshooting

### Test Fails to Find Headers

Ensure the test support library is linked:

```cmake
target_link_libraries(my_test PRIVATE emulator_test_support)
```

### MPI Tests Hang

Check that all ranks participate in collective operations and that
`MPI_Barrier` is used before/after critical sections.

### Tests Pass Locally but Fail in CI

- Check for hardcoded paths
- Verify random seeds are set for reproducibility
- Check for timing-dependent assertions
