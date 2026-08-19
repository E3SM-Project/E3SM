#!/bin/bash
# Profile one standalone P3 SP run on a Perlmutter GPU node.
#
# This intentionally does not run the normal timing or correctness matrix.
# Use PROFILE_TOOL=nsys for a low-overhead timeline, or PROFILE_TOOL=ncu for
# kernel metrics. The executable must already have been built on the login
# node with ./build-p3-profile-gpu-sp.sh. The profiling build is kept in a
# separate directory so it cannot overwrite the normal testing executable.
#
# Examples:
#   sbatch sbatch-p3-profile-pm-gpu-sp.sh
#   PROFILE_TOOL=ncu NCU_KERNEL_NAME='regex:.*p3_.*_disp.*' \
#     NCU_LAUNCH_SKIP=0 NCU_LAUNCH_COUNT=0 \
#     sbatch sbatch-p3-profile-pm-gpu-sp.sh

#SBATCH -A e3sm
#SBATCH -q debug
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH --gpus=1
#SBATCH -t 00:30:00
#SBATCH -J p3_profile_pm-gpu-sp
#SBATCH -o p3-profile-results/o.p3_profile.%j.txt

set -euo pipefail

if [ -n "${NERSC_HOST:-}" ] && [ "$NERSC_HOST" != perlmutter ]; then
  echo "Error: NERSC_HOST=${NERSC_HOST}; expected perlmutter." >&2
  exit 1
fi

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
TEST_EXE="${REPO_ROOT}/build_p3_profile_pm-gpu-SP/src/physics/p3/tests/p3_run_and_cmp"
PROFILE_TOOL="${PROFILE_TOOL:-nsys}"
SCRATCH_ROOT="${SCRATCH:-/pscratch/sd/${USER:0:1}/${USER}}"
RESULTS_DIR="${RESULTS_DIR:-${SCRATCH_ROOT}/e3sm/p3-profile-results/pm-gpu-sp/${SLURM_JOB_ID:-manual}}"
BASELINE_DIR="${RESULTS_DIR}/baseline"
NCOL="${NCOL:-20000}"
NLEV="${NLEV:-128}"
NSTEPS="${NSTEPS:-1}"
WARMUP="${WARMUP:-1}"
DT="${DT:-300}"
PREDICT_NC="${PREDICT_NC:-yes}"
PRESCRIBED_CCN="${PRESCRIBED_CCN:-yes}"
PROFILE_REPETITIONS="${PROFILE_REPETITIONS:-1}"
NCU_LAUNCH_SKIP="${NCU_LAUNCH_SKIP:-0}"
# Zero means unrestricted collection; do not pass --launch-count to ncu.
NCU_LAUNCH_COUNT="${NCU_LAUNCH_COUNT:-0}"
# Empty means profile all eligible kernels; set NCU_KERNEL_NAME to filter.
NCU_KERNEL_NAME="${NCU_KERNEL_NAME-}"

if [ ! -x "$TEST_EXE" ]; then
  echo "Error: missing executable: ${TEST_EXE}" >&2
  echo "Build it first with ./build-p3-profile-gpu-sp.sh" >&2
  exit 1
fi
if ! command -v "$PROFILE_TOOL" >/dev/null 2>&1; then
  echo "Error: ${PROFILE_TOOL} is not available in this environment." >&2
  exit 1
fi

mkdir -p "$SCRATCH_ROOT"
if [ ! -w "$SCRATCH_ROOT" ]; then
  echo "Error: scratch directory is not writable: ${SCRATCH_ROOT}" >&2
  echo "Set RESULTS_DIR to a writable scratch location." >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR" "$BASELINE_DIR"
cd "$RESULTS_DIR"
export TMPDIR="${PROFILE_TMPDIR:-${RESULTS_DIR}/tmp}"
mkdir -p "$TMPDIR"
ARGS=(-g -b "$BASELINE_DIR" -s "$NSTEPS" -w "$WARMUP" -dt "$DT" -i "$NCOL" -k "$NLEV"
      --predict-nc "$PREDICT_NC" --prescribed-ccn "$PRESCRIBED_CCN"
      -r "$PROFILE_REPETITIONS")

echo "Profiling ${TEST_EXE} with ${PROFILE_TOOL}"
echo "Results: ${RESULTS_DIR}"
echo "Arguments: ${ARGS[*]}"
echo "NCU kernel filter: ${NCU_KERNEL_NAME:-<none>}"
echo "CUDA toolkit/runtime: ${CUDA_HOME:-unknown}"
command -v nsys >/dev/null 2>&1 && nsys --version || true
nvidia-smi --query-gpu=driver_version --format=csv,noheader || true

DCGM_PAUSED=0
DCGM_REPAUSE_PID=""
pause_dcgm_for_ncu() {
  if ! command -v dcgmi >/dev/null 2>&1; then
    echo "WARNING: dcgmi not found; ncu may fail because DCGM owns the profiling counters." >&2
    return
  fi
  echo "Pausing DCGM profiling counters for ncu..."
  if dcgmi profile --pause >/dev/null 2>&1; then
    DCGM_PAUSED=1
    (
      while sleep 120; do
        dcgmi profile --pause >/dev/null 2>&1 || true
      done
    ) &
    DCGM_REPAUSE_PID=$!
  else
    echo "WARNING: dcgmi profile --pause failed; ncu may not get exclusive counter access." >&2
  fi
}
resume_dcgm() {
  if [ -n "$DCGM_REPAUSE_PID" ]; then
    kill "$DCGM_REPAUSE_PID" >/dev/null 2>&1 || true
    wait "$DCGM_REPAUSE_PID" 2>/dev/null || true
  fi
  if [ "$DCGM_PAUSED" -eq 1 ]; then
    dcgmi profile --resume >/dev/null 2>&1 || true
  fi
}

case "$PROFILE_TOOL" in
  nsys)
    nsys profile --force-overwrite=true --sample=none --cpuctxsw=none \
      --trace=cuda,nvtx,osrt --stats=true -o "${RESULTS_DIR}/p3_profile" \
      "$TEST_EXE" "${ARGS[@]}" 2>&1 | tee nsys.log
    echo "Report: ${RESULTS_DIR}/p3_profile.nsys-rep"
    ;;
  ncu)
    pause_dcgm_for_ncu
    trap resume_dcgm EXIT
    # NCU may replay each selected launch internally; launch-count limits
    # which matching launches are collected.
    NCU_ARGS=(--set full --target-processes all --kernel-name-base demangled
              --launch-skip "$NCU_LAUNCH_SKIP"
              -o p3_profile)
    if [ "$NCU_LAUNCH_COUNT" -gt 0 ]; then
      NCU_ARGS+=(--launch-count "$NCU_LAUNCH_COUNT")
    fi
    if [ -n "$NCU_KERNEL_NAME" ]; then
      NCU_ARGS+=(--kernel-name "$NCU_KERNEL_NAME")
    fi
    ncu "${NCU_ARGS[@]}" "$TEST_EXE" "${ARGS[@]}" 2>&1 | tee ncu.log
    echo "Report: ${RESULTS_DIR}/p3_profile.ncu-rep"
    ;;
  *)
    echo "Error: PROFILE_TOOL must be nsys or ncu, got ${PROFILE_TOOL}" >&2
    exit 1
    ;;
esac
