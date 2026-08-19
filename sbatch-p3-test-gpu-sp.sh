#!/bin/bash
#SBATCH -A e3sm
#SBATCH -q debug
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH --gpus=4
#SBATCH -t 00:30:00
#SBATCH -J p3-test-gpu-sp
#SBATCH -o p3-test-results/o.p3-test-gpu-sp.%j.txt
#
# Standalone P3 GPU test in single precision.  The machine is selected from
# NERSC_HOST so the same script can be submitted on Perlmutter or Alvarez.
# Build first on the login node with ./build-p3-test-gpu-sp.sh.
set -euo pipefail

case "${NERSC_HOST:-}" in
  perlmutter) MACHINE_NAME=pm-gpu ;;
  alvarez) MACHINE_NAME=alvarez-gpu ;;
  *)
    echo "Error: NERSC_HOST must be perlmutter or alvarez; got '${NERSC_HOST:-unset}'." >&2
    exit 1
    ;;
esac
export MACHINE_NAME
export NPAR="${NPAR:-1}"
export NCOL="${NCOL:-20000}"
export NLEV="${NLEV:-128}"
export NSTEPS="${NSTEPS:-6}"
export DT="${DT:-300}"
export WARMUP="${WARMUP:-3}"
export REPEAT="${REPEAT:-20}"
export VERIFY_MODE="${VERIFY_MODE:-none}"
export NRUNS="${NRUNS:-2}"
export TOL="${TOL:-0}"

if [ "$MACHINE_NAME" = alvarez-gpu ]; then
  CUDA_LIB_DIR="${CUDA_LIB_DIR:-/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda/12.9/lib64}"
  export LD_LIBRARY_PATH="${CUDA_LIB_DIR}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
exec "${REPO_ROOT}/run_p3_test_${MACHINE_NAME}-sp.sh" "$MACHINE_NAME"
