#!/bin/bash
# Shared build helper for the standalone p3_run_and_cmp EAMxx test.
#
# Configures and builds (login-node friendly: no GPU/many-core resources
# needed) the p3_run_and_cmp executable for a given machine + precision.
# Intended to be called by the thin build_p3_test_*.sh wrappers; keeping the
# build step separate from run_p3_test.sh means you can build once on the
# login node and then run (possibly many times, possibly in parallel) on a
# compute node/allocation without ever touching compute-node resources for
# compilation.
#
# Usage:
#   build_p3_test.sh <machine-name> <SP|DP> <build-dir> [Release|Debug] [ON|OFF]
#
# All standalone P3 builds use the small-kernel dispatch path.  This keeps the
# executable used for testing and profiling consistent across machines.
#
set -e

MACHINE_NAME="$1"
PRECISION="$2"   # SP or DP
BUILD_DIR="$3"
BUILD_TYPE="${4:-Release}"
CUDA_PROFILER="${5:-OFF}"

if [ -z "$MACHINE_NAME" ] || [ -z "$PRECISION" ] || [ -z "$BUILD_DIR" ]; then
  echo "Usage: build_p3_test.sh <machine-name> <SP|DP> <build-dir> [Release|Debug] [ON|OFF]" >&2
  exit 1
fi

case "$PRECISION" in
  SP) DOUBLE_PRECISION=OFF ;;
  DP) DOUBLE_PRECISION=ON ;;
  *) echo "Error: precision must be SP or DP, got '${PRECISION}'" >&2; exit 1 ;;
esac

case "$CUDA_PROFILER" in
  ON|OFF) ;;
  *) echo "Error: CUDA profiler option must be ON or OFF, got '${CUDA_PROFILER}'" >&2; exit 1 ;;
esac

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
EAMXX_DIR="${REPO_ROOT}/components/eamxx"

if [ ! -d "$EAMXX_DIR" ]; then
    echo "Error: could not locate components/eamxx under repo root (${REPO_ROOT})." >&2
    exit 1
fi

# CIME now requires Python 3.9 or newer.  On Alvarez, an inherited Anaconda
# path can place Python 3.8 ahead of the loaded cray-python module, so make
# the module's interpreter selection explicit before invoking CIME.
if [[ "$MACHINE_NAME" == alvarez-* ]] && [ -x /opt/cray/pe/python/3.12.12/bin/python3 ]; then
    export PATH="/opt/cray/pe/python/3.12.12/bin:${PATH}"
    hash -r 2>/dev/null || true
fi

echo "=========================================================="
echo " Loading modules for ${MACHINE_NAME}"
echo "=========================================================="
export SHELL=/bin/bash
eval $(${EAMXX_DIR}/scripts/eamxx-env-cmd ${MACHINE_NAME})

echo "=========================================================="
echo " Setting up build directory for EAMxx P3 standalone test"
echo " (${MACHINE_NAME}, ${PRECISION}) -> ${BUILD_DIR}"
echo "=========================================================="

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# 1. Configure via CMake
# We use the machine-specific initialization cache file to get compiler/MPI/Kokkos settings right.
# We also set options to mimic a typical ne256 run (Release build, enable P3, etc.)
# Compiler names are machine-specific. For example, Cray systems like
# pm-cpu/pm-gpu use CC/cc/ftn, while module-based systems like vista-gh use
# mpicxx/mpicc/mpifort. Query the EAMxx machine spec instead of hardcoding.
CXX_COMPILER="$("${EAMXX_DIR}/scripts/query-eamxx" "${MACHINE_NAME}" cxx_compiler)"
C_COMPILER="$("${EAMXX_DIR}/scripts/query-eamxx" "${MACHINE_NAME}" c_compiler)"
FTN_COMPILER="$("${EAMXX_DIR}/scripts/query-eamxx" "${MACHINE_NAME}" f90_compiler)"

echo "Running CMake..."
cmake -DCMAKE_CXX_COMPILER="${CXX_COMPILER}" -DCMAKE_C_COMPILER="${C_COMPILER}" -DCMAKE_Fortran_COMPILER="${FTN_COMPILER}" -C "${EAMXX_DIR}/cmake/machine-files/${MACHINE_NAME}.cmake" \
      -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
      ${BUILD_TYPE:+-DKokkos_ENABLE_DEBUG=\"$([ "${BUILD_TYPE}" = Debug ] && echo ON || echo OFF)\"} \
      ${BUILD_TYPE:+-DKokkos_ENABLE_DEBUG_BOUNDS_CHECK=\"$([ "${BUILD_TYPE}" = Debug ] && echo ON || echo OFF)\"} \
      -DSCREAM_DOUBLE_PRECISION=${DOUBLE_PRECISION} \
      -DSCREAM_SMALL_KERNELS=ON \
      -DSCREAM_P3_SMALL_KERNELS=ON \
      -DSCREAM_ENABLE_CUDA_PROFILER=${CUDA_PROFILER} \
      -DSCREAM_TEST_MAX_THREADS=1 \
      -DSCREAM_ENABLE_BASELINE_TESTS=OFF \
      ${NETCDF_PATH:+-DNetCDF_PATH="${NETCDF_PATH}"} \
      ${PNETCDF_PATH:+-DPnetCDF_PATH="${PNETCDF_PATH}"} \
      "${EAMXX_DIR}"

# 2. Build the specific test executable
# p3_run_and_cmp is typically the heaviest and closest to an actual run loop for P3 tests.
echo "=========================================================="
echo " Building p3_run_and_cmp test..."
echo "=========================================================="
make -j8 p3_run_and_cmp

EXE="${BUILD_DIR}/src/physics/p3/tests/p3_run_and_cmp"
if [ ! -x "$EXE" ]; then
  echo "Error: expected executable not found at ${EXE}" >&2
  exit 1
fi

echo "=========================================================="
echo " Build complete: ${EXE}"
echo " Small P3 kernels: ON"
echo " (build performed on login node; use run_p3_test_*.sh on a"
echo "  compute node/allocation to actually execute the test)"
echo "=========================================================="
