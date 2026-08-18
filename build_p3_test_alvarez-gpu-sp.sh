#!/bin/bash
# Build (login-node only, no compute resources needed) the p3_run_and_cmp
# standalone test for the Alvarez GPU machine, SP.
#
# See scripts_p3_test/build_p3_test.sh for details. After this completes,
# use run_p3_test_alvarez-gpu-sp.sh on a compute node/allocation to execute it.
set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACHINE_NAME="alvarez-gpu"
BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-SP"
"${REPO_ROOT}/scripts_p3_test/build_p3_test.sh" "$MACHINE_NAME" SP "$BUILD_DIR"
