#!/bin/bash
# Run the Alvarez GPU standalone P3 test in single precision.
set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACHINE_NAME="${MACHINE_NAME:-alvarez-gpu}"
BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-SP"
TEST_EXE="${BUILD_DIR}/src/physics/p3/tests/p3_run_and_cmp"

if [ ! -x "$TEST_EXE" ]; then
  echo "Error: SP executable is missing or not executable: ${TEST_EXE}" >&2
  echo "Build it first with ./build_p3_test_alvarez-gpu-sp.sh" >&2
  exit 1
fi

if [[ "$MACHINE_NAME" == alvarez-* ]]; then
  CUDA_LIB_DIR="${CUDA_LIB_DIR:-/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda/12.9/lib64}"
  export LD_LIBRARY_PATH="${CUDA_LIB_DIR}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi

exec "${REPO_ROOT}/scripts_p3_test/run_p3_test.sh" "$TEST_EXE" "${MACHINE_NAME}-sp"
