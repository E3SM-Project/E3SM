#!/bin/bash
# Build (login-node only, no compute resources needed) the p3_run_and_cmp
# standalone test for the Alvarez GPU machine, DP.
#
# See scripts_p3_test/build_p3_test.sh for details. After this completes,
# use run_p3_test_alvarez-gpu-dp.sh on a compute node/allocation to execute it.
set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACHINE_NAME="alvarez-gpu"
BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-DP"
"${REPO_ROOT}/scripts_p3_test/build_p3_test.sh" "$MACHINE_NAME" DP "$BUILD_DIR"
