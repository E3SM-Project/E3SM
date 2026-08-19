#!/bin/bash
#SBATCH -A e3sm
#SBATCH -q regular
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH --gpus=4
#SBATCH -t 02:00:00
#SBATCH -J sbatch_p3_paired11NCOL_alvarez-gpu
#SBATCH -o alvarez-gpu-p3_test_results/o.sbatch_p3_paired11NCOL_alvarez-gpu.%j.txt
#
# Paired SP/DP NCOL timing sweep for pm-gpu.  Each NCOL point runs SP and DP
# sequentially on one GPU, so the two precision measurements share the same
# GPU and do not run concurrently with each other.  Four different points may
# run concurrently, one per GPU.
#
# Build both executables on the login node before submitting:
#   NERSC_HOST=perlmutter ./build-p3-test-gpu-sp.sh
#   ./build-p3-test-gpu-dp.sh
#
# This legacy script intentionally keeps its historical underscore name.
# Submit from the repository root:
#   sbatch sbatch_p3_paired11NCOL_alvarez-gpu.sh
set -e
MACHINE_NAME="${MACHINE_NAME:-alvarez-gpu}"

# Alvarez does not provide the Environment Modules command.  The executables
# were built with the CUDA 12.9 toolkit, so use its runtime library explicitly.
CUDA_LIB_DIR="${CUDA_LIB_DIR:-/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda/12.9/lib64}"
export LD_LIBRARY_PATH="${CUDA_LIB_DIR}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

# Match the prior SP sweep, including its logarithmic tail to 40000.  Add
# targeted points in the low-NCOL plateau and the noisier transition region.
# The two passes below make this a two-pass sweep.
# Each completed SP/DP point also writes total elapsed time and average
# measured kernel time to the Slurm output file for live monitoring.
NCOL_LIST=(
  100 150 200 250 300 350
  400 450 500 550 600 650 700 800 850 900 1000
  1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
  2342 2741 3210 3758 4400 5151 6031 7060 8266 9678 11331 13266 15531
  18184 21289 24925 29182 34165 40000
)

export NLEV=128
export NSTEPS=6
export DT=300
export WARMUP=1
export REPEAT=5
export VERIFY_MODE=none
export NRUNS=0
export TOL=0
export PROFILE=off
export NPAR=1
export SWEEP_PASSES=1
export GPU_COUNT=4

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SP_EXE="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-SP/src/physics/p3/tests/p3_run_and_cmp"
DP_EXE="${REPO_ROOT}/build_p3_test_${MACHINE_NAME}-DP/src/physics/p3/tests/p3_run_and_cmp"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

check_executable() {
  local label="$1" exe="$2"
  [ -x "$exe" ] || die "${label} executable is missing or not executable: ${exe}"
  [ -s "$exe" ] || die "${label} executable is empty: ${exe}"
  file "$exe" | grep -q 'ELF 64-bit.*x86-64' || \
    die "${label} executable does not look like an x86-64 ELF binary: ${exe}"
}

check_executable SP "$SP_EXE"
check_executable DP "$DP_EXE"

echo "CUDA_LIB_DIR=${CUDA_LIB_DIR}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
ldd "$SP_EXE" | grep -E 'cuda|not found' || true

# Both builds must have the same basic binary and runtime shape.  Precision
# changes are expected to alter the executable hash, but a differing ABI,
# architecture, or shared-library set usually means one build is stale or was
# made for a different target.
SP_FILE_INFO="$(file "$SP_EXE")"
DP_FILE_INFO="$(file "$DP_EXE")"
SP_ELF_INFO="$(readelf -h "$SP_EXE" | awk -F: '/Class:|Data:|Machine:|OS\/ABI:/ {gsub(/^ +/,"",$2); print $1 "=" $2}')"
DP_ELF_INFO="$(readelf -h "$DP_EXE" | awk -F: '/Class:|Data:|Machine:|OS\/ABI:/ {gsub(/^ +/,"",$2); print $1 "=" $2}')"
[ "$SP_ELF_INFO" = "$DP_ELF_INFO" ] || die "SP and DP ELF headers are inconsistent"
# Record dependency output for diagnostics, but do not require exact matches.
# SP and DP can resolve the same libraries through different loader paths or
# report environment-dependent details such as CUDA runtime versions.
SP_LDD="$(ldd "$SP_EXE" 2>&1 | sed -E 's/ \(0x[0-9a-fA-F]+\)//g')"
DP_LDD="$(ldd "$DP_EXE" 2>&1 | sed -E 's/ \(0x[0-9a-fA-F]+\)//g')"
echo "SP shared-library dependencies:"
echo "$SP_LDD"
echo "DP shared-library dependencies:"
echo "$DP_LDD"

RESULTS_ROOT="${REPO_ROOT}/alvarez-gpu-p3_test_results"
SWEEP_ROOT="${RESULTS_ROOT}/ncol_sweep_pm-gpu-sp-dp-paired"
mkdir -p "$SWEEP_ROOT"
JOB_ID="${SLURM_JOB_ID:-local}"
METADATA="${SWEEP_ROOT}/ncol_sweep_pm-gpu-sp-dp-paired.${JOB_ID}.metadata.txt"
{
  echo "timestamp=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  echo "hostname=$(hostname)"
  echo "slurm_job_id=${JOB_ID}"
  echo "slurm_job_name=${SLURM_JOB_NAME:-local}"
  echo "nlev=${NLEV}"
  echo "nsteps=${NSTEPS}"
  echo "warmup=${WARMUP}"
  echo "repeat=${REPEAT}"
  echo "passes=${SWEEP_PASSES}"
  echo "gpu_count=${GPU_COUNT}"
  echo "sp_exe=${SP_EXE}"
  echo "sp_sha256=$(sha256sum "$SP_EXE" | awk '{print $1}')"
  echo "dp_exe=${DP_EXE}"
  echo "dp_sha256=$(sha256sum "$DP_EXE" | awk '{print $1}')"
  echo "sp_file=${SP_FILE_INFO}"
  echo "dp_file=${DP_FILE_INFO}"
  echo "nvidia_smi_query:"
  nvidia-smi --query-gpu=index,name,memory.total,clocks.current.graphics,pstate --format=csv 2>&1 || true
  echo "nvidia_smi_list:"
  nvidia-smi -L 2>&1 || true
} > "$METADATA"
cat "$METADATA"

run_precision() {
  local precision="$1" ncol="$2" pass="$3" gpu="$4"
  local exe tag run_root log start_time end_time total_sec run_dir kernel_avg
  if [ "$precision" = DP ]; then exe="$DP_EXE"; else exe="$SP_EXE"; fi
  tag="${MACHINE_NAME}-${precision,,}-ncol${ncol}-pass${pass}-gpu${gpu}"
  run_root="${SWEEP_ROOT}/${tag}"
  log="${SWEEP_ROOT}/o.${tag}.${JOB_ID}.txt"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] start ${precision} NCOL=${ncol} pass=${pass} GPU=${gpu}"
  start_time=$(date +%s.%N)
  set +e
  NCOL="$ncol" CUDA_VISIBLE_DEVICES="$gpu" HIP_VISIBLE_DEVICES="$gpu" \
    RESULTS_DIR="$run_root" \
    "$REPO_ROOT/scripts_p3_test/run_p3_test.sh" "$exe" "$tag" \
    > "$log" 2>&1
  local rc=$?
  set -e
  end_time=$(date +%s.%N)
  total_sec=$(awk -v a="$start_time" -v b="$end_time" 'BEGIN {printf "%.3f", b-a}')
  run_dir=$(find "$run_root" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort | tail -1)
  kernel_avg="NA"
  if [ -n "$run_dir" ] && [ -f "${run_dir}/timing.log" ]; then
    kernel_avg=$(awk '
      /^Measured steps: mean=/ {
        split($3, a, "=")
        sum += a[2]
        count++
      }
      END {
        if (count > 0) printf "%.6e", sum / count
        else print "NA"
      }
    ' "${run_dir}/timing.log")
  fi
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] done ${precision} NCOL=${ncol} pass=${pass} GPU=${gpu} status=$([ "$rc" -eq 0 ] && echo PASS || echo FAIL) total_sec=${total_sec} avg_kernel_sec_per_step=${kernel_avg}"
  return "$rc"
}

run_pair() {
  local ncol="$1" pass="$2" gpu="$3" first second rc=0
  if [ "$pass" -eq 1 ]; then
    first=DP
    second=SP
  else
    first=SP
    second=DP
  fi
  run_precision "$first" "$ncol" "$pass" "$gpu" || rc=1
  run_precision "$second" "$ncol" "$pass" "$gpu" || rc=1
  return "$rc"
}

# Pass 1 uses DP-then-SP and pass 2 uses SP-then-DP.  Rotating the GPU in the
# second pass separates precision/order effects from a fixed GPU identity.
FAIL=0
for pass in $(seq 1 "$SWEEP_PASSES"); do
  PIDS=()
  point_index=0
  for ncol in "${NCOL_LIST[@]}"; do
    gpu=$(( (point_index + pass - 1) % GPU_COUNT ))
    run_pair "$ncol" "$pass" "$gpu" &
    PIDS+=("$!")
    point_index=$((point_index + 1))
    if [ "${#PIDS[@]}" -ge "$GPU_COUNT" ]; then
      for pid in "${PIDS[@]}"; do wait "$pid" || FAIL=1; done
      PIDS=()
    fi
  done
  for pid in "${PIDS[@]}"; do wait "$pid" || FAIL=1; done
done

CSV="${SWEEP_ROOT}/ncol_sweep_pm-gpu-sp-dp-paired.${JOB_ID}.csv"
echo "precision,ncol,pass,gpu,case,kernel_time_sec_per_step" > "$CSV"
for pass in $(seq 1 "$SWEEP_PASSES"); do
  point_index=0
  for ncol in "${NCOL_LIST[@]}"; do
    gpu=$(( (point_index + pass - 1) % GPU_COUNT ))
    point_index=$((point_index + 1))
    for precision in dp sp; do
      tag="${MACHINE_NAME}-${precision}-ncol${ncol}-pass${pass}-gpu${gpu}"
      run_dir=$(find "${SWEEP_ROOT}/${tag}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort | tail -1)
      case_no=0
      if [ -n "$run_dir" ] && [ -f "${run_dir}/timing.log" ]; then
        grep -E '^Measured steps: mean=' "${run_dir}/timing.log" | while read -r line; do
          value=$(echo "$line" | awk '{split($3,a,"="); print a[2]}')
          echo "${precision^^},${ncol},${pass},${gpu},${case_no},${value}" >> "$CSV"
          case_no=$((case_no + 1))
        done
      fi
    done
  done
done

echo "Paired SP/DP sweep done. CSV: ${CSV}"
cat "$CSV"
if [ "$FAIL" -ne 0 ]; then
  die "one or more paired runs failed; inspect per-run logs under ${SWEEP_ROOT}"
fi
