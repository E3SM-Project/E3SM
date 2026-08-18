#!/bin/bash
# Build the standalone p3_run_and_cmp test for Perlmutter GPU nodes, DP.
# Run this on a Perlmutter login node, then use the matching sbatch script.
set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACHINE_NAME="pm-gpu"

if [ -n "${NERSC_HOST:-}" ] && [ "$NERSC_HOST" != perlmutter ]; then
  echo "Error: NERSC_HOST=${NERSC_HOST}; expected perlmutter for pm-gpu." >&2
  exit 1
fi

BUILD_DIR="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-DP"
"${REPO_ROOT}/scripts_p3_test/build_p3_test.sh" "$MACHINE_NAME" DP "$BUILD_DIR"
