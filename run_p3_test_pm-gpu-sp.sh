#!/bin/bash
# Run the Perlmutter GPU standalone P3 test in single precision.
set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACHINE_NAME="pm-gpu"
BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-SP"
TEST_EXE="${BUILD_DIR}/src/physics/p3/tests/p3_run_and_cmp"

if [ ! -x "$TEST_EXE" ]; then
  echo "Error: SP executable is missing or not executable: ${TEST_EXE}" >&2
  echo "Build it first with ./build_p3_test_pm-gpu-sp.sh" >&2
  exit 1
fi

exec "${REPO_ROOT}/scripts_p3_test/run_p3_test.sh" "$TEST_EXE" "${MACHINE_NAME}-sp"
