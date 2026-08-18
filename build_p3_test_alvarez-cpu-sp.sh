#!/bin/bash
# Build the standalone P3 test for the Alvarez CPU partition, SP.
set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACHINE_NAME="${MACHINE_NAME:-alvarez-cpu}"
BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-SP"
"${REPO_ROOT}/scripts_p3_test/build_p3_test.sh" "$MACHINE_NAME" SP "$BUILD_DIR"
