#!/bin/bash
#SBATCH -A e3sm
#SBATCH -q debug
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH --gpus=4
#SBATCH -t 00:30:00
#SBATCH -J p3_test_pm-gpu-sp
#SBATCH -o pm-gpu-p3_test_results/o.p3_test_pm-gpu-sp.%j.txt
#
# Perlmutter GPU batch job for the standalone P3 test in SP.
# Build first on the login node with ./build_p3_test_pm-gpu-sp.sh.
set -e

if [ -n "${NERSC_HOST:-}" ] && [ "$NERSC_HOST" != perlmutter ]; then
  echo "Error: NERSC_HOST=${NERSC_HOST}; expected perlmutter for pm-gpu." >&2
  exit 1
fi

export MACHINE_NAME=pm-gpu
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

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
exec "${REPO_ROOT}/run_p3_test_pm-gpu-sp.sh"
