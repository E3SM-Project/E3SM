#!/bin/bash
# Build a debuggable SP P3 executable on an Alvarez login node.
set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACHINE_NAME="alvarez-gpu"
BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-SP-debug"
"${REPO_ROOT}/scripts_p3_test/build_p3_test.sh" "$MACHINE_NAME" SP "$BUILD_DIR" Debug
