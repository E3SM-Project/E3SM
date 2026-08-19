#!/bin/bash
# Build the standalone p3_run_and_cmp test in double precision.
# Run this on the login node for the NERSC system being used.
set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "${NERSC_HOST:-}" in
  perlmutter)
    MACHINE_NAME="pm-gpu"
    BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-DP"
    ;;
  alvarez)
    MACHINE_NAME="alvarez-gpu"
    BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-DP"
    ;;
  *)
    echo "Error: NERSC_HOST must be perlmutter or alvarez; got '${NERSC_HOST:-unset}'." >&2
    exit 1
    ;;
esac

"${REPO_ROOT}/scripts_p3_test/build_p3_test.sh" \
  "$MACHINE_NAME" DP "$BUILD_DIR"
