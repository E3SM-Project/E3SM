# Standalone Build Guide

Build emulator_comps independently from CIME for development and testing.

## Quick Start

Use the `test` script for an automated build experience:

```console
cd components/emulator_comps

# Build with stub backend (default, no dependencies)
./test rebuild

# Build with LibTorch backend (auto-detects from pytorch module)
./test rebuild --backend libtorch

# Run tests only
./test test
```

## Test Script Usage

```console
./test [action] [options]

Actions:
  build       Configure (if needed) and build (default)
  rebuild     Clean, configure, and build
  test        Run tests only (must be built first)
  configure   CMake configuration only
  clean       Remove build directory

Options:
  -h, --help              Show help
  -b, --backend BACKEND   Inference backend: stub (default), libtorch
  -j, --jobs N            Parallel jobs (default: 16)
  -jN                     Shorthand (e.g., -j8)
  -t, --test-level N      Test level (1=unit, 2=component, 3=all)
  -q, --quiet             Minimal output (logs to files)
  -v, --verbose           Full output
  --no-tests              Skip tests
  --mpi                   Enable MPI tests (requires compute node)
  --debug                 Debug build
```

### Backends

| Backend | Description | Dependencies |
| ------- | ----------- | ------------ |
| `stub` | No-op backend for testing | None |
| `libtorch` | LibTorch C++ inference | pytorch module (Perlmutter) |

### Examples

```console
# Build with stub backend (no dependencies)
./test rebuild --backend stub

# Build with LibTorch on Perlmutter
./test rebuild --backend libtorch

# Quick build with 8 cores, quiet mode
./test rebuild -j8 -q

# Unit tests only
./test rebuild -t1

# Full test suite on compute node (with MPI)
salloc -N 1 -t 30:00 -C cpu
./test rebuild --mpi
```

> **Note:** MPI tests require `srun` which only works on compute nodes.
> The script skips MPI tests by default for login node compatibility.

## Script Structure

The test script is modular with machine-specific setup in `scripts/`:

```console
components/emulator_comps/
├── test                    # Main entry point
└── scripts/
    ├── common.sh           # Shared functions (logging, usage)
    ├── perlmutter.sh       # NERSC Perlmutter setup + LibTorch detection
    └── generic.sh          # Generic Linux setup
```

## Manual Build

### On Perlmutter (NERSC)

```console
# Load required modules (ORDER MATTERS: HDF5 before NetCDF)
module load PrgEnv-gnu
module load craype-x86-milan
module load cmake
module load cray-hdf5
module load cray-netcdf
module load cray-parallel-netcdf

# Optional: Load pytorch for LibTorch backend
module load pytorch/2.6.0
export Torch_DIR=$(python -c "import torch; print(torch.__path__[0])")/share/cmake/Torch

# Navigate to emulator_comps
cd $SCRATCH/e3sm/components/emulator_comps

# Remove any previous build
rm -rf build && mkdir build && cd build

# Configure standalone build with tests
cmake .. \
  -DEMULATOR_STANDALONE_BUILD=ON \
  -DEMULATOR_BUILD_TESTS=ON \
  -DEATM_INFERENCE_BACKEND=stub \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_Fortran_COMPILER=ftn

# Build
make -j16

# Run tests
ctest --output-on-failure
```

### On Generic Linux

```console
# Ensure MPI compilers are in PATH
export PATH=/path/to/mpi/bin:$PATH

cd components/emulator_comps
mkdir build && cd build

cmake .. \
  -DEMULATOR_STANDALONE_BUILD=ON \
  -DEMULATOR_BUILD_TESTS=ON \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90

make -j$(nproc)
ctest --output-on-failure
```

## Requirements

**Compilers:**

- C compiler (gcc 8+, Intel 19+, or Clang 10+)
- C++ compiler with C++17 support
- Fortran compiler (gfortran 8+ or ifort 19+)

**Libraries:**

- MPI (OpenMPI, MPICH, Intel MPI, or Cray MPICH)
- NetCDF-C (for SCORPIO)

**Optional:**

- PyTorch/LibTorch (for libtorch backend)

**Tools:**

- CMake 3.16+
- Git (for submodules)

## What Gets Built

In standalone mode, CMake automatically builds dependencies from E3SM externals:

| Dependency | Source Location | Purpose |
| ---------- | --------------- | ------- |
| SCORPIO | `externals/scorpio` | Parallel I/O |
| yaml-cpp | `externals/ekat/extern/yaml-cpp` | YAML config parsing |
| Catch2 | `externals/ekat/extern/Catch2` | Test framework |

## Build Options

| Option | Default | Description |
| ------ | ------- | ----------- |
| `EMULATOR_STANDALONE_BUILD` | OFF | Enable standalone build mode |
| `EMULATOR_BUILD_TESTS` | OFF | Build test executables |
| `EATM_INFERENCE_BACKEND` | stub | Inference backend (stub, libtorch) |
| `EMULATOR_TEST_LEVEL` | 3 | Test level (1=unit, 2=component, 3=system) |
| `EMULATOR_SKIP_MPI_TESTS` | OFF | Skip MPI tests (for login nodes) |
| `EMULATOR_HAS_LIBTORCH` | OFF | Enable LibTorch (auto-set by backend) |

## Directory Structure After Build

```console
build/
├── logs/                 # Build and test logs
│   ├── cmake_configure.log
│   ├── build.log
│   ├── test.log
│   └── build_summary.txt
├── externals/
│   ├── scorpio/          # Built SCORPIO library
│   └── yaml-cpp/         # Built yaml-cpp library
├── common/
│   ├── libemulator_common.a
│   └── tests/
├── eatm/
│   ├── libeatm.a
│   └── tests/
└── CTestTestfile.cmake
```

## Troubleshooting

### SCORPIO Build Fails

**Symptom:** Errors about NetCDF not found

**Solution:** Ensure NetCDF is available:

```console
nc-config --version
nc-config --includedir
cmake .. -DNetCDF_C_PATH=/path/to/netcdf
```

### LibTorch Not Found

**Symptom:** `Torch_DIR not set` error

**Solution on Perlmutter:**

```console
module load pytorch/2.6.0
```

**Solution on generic Linux:**

```console
export Torch_DIR=/path/to/libtorch/share/cmake/Torch
```

### MPI Tests Fail with "srun error"

**Symptom:** `srun: error: No architecture specified`

**Solution:** MPI tests require a compute node:

```console
# Option 1: Skip MPI tests (default)
./test rebuild

# Option 2: Run on compute node
salloc -N 1 -t 30:00 -C cpu
./test rebuild --mpi
```

### yaml-cpp Not Found

**Solution:** Initialize EKAT submodule:

```console
cd $E3SM_ROOT
git submodule update --init --recursive externals/ekat
```
