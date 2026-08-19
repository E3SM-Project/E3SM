#!/bin/bash
# Shared runner for the standalone p3_run_and_cmp EAMxx test.
#
# This is called by the run_p3_test_*.sh wrappers for supported machines.
# against an *already-built* executable (see build_p3_test_*.sh /
# build_p3_test.sh, which are meant to be run on the login node) -- it is the
# one shared runner, so any knob/behavior added here applies uniformly to all
# 4. It is meant to be invoked on a compute node/allocation (interactive
# salloc or inside an sbatch job), since it's the part that actually burns
# compute-node resources.
#
# It runs the test with configurable "difficulty" knobs (problem size,
# timestep count/length, repeat count for timing), saves timing results (both
# the P3 kernel's own internal std::chrono::steady_clock timing, see REPEAT
# below, and wall-clock time of every exe invocation), and checks full-field
# hashes so that separate runs (same machine, different machine, before/after
# a code change, etc.) can be checked for bit-for-bit reproducibility without
# keeping large baseline files by default.
#
# Usage:
#   run_p3_test.sh <path-to-p3_run_and_cmp-exe> <tag>
#
# <tag> is a short label (e.g. pm-gpu-dp) used to namespace results/baselines
# on disk so that runs from different build variants don't clobber each
# other. It is also used to decide GPU vs CPU pinning for NPAR (see below):
# Tags containing gpu (and vista-gh) are treated as GPU, everything else as CPU.
#
# All knobs are controlled via environment variables (all optional):
#   NCOL           number of columns              (default 3)
#   NLEV           number of vertical levels       (default 72)
#   NSTEPS         number of timesteps             (default 6)
#   WARMUP         warmup timesteps before measured timing steps (default 3)
#   DT             timestep length in seconds      (default 300)
#   PREDICT_NC     yes|no|both                     (default both)
#   PRESCRIBED_CCN yes|no|both                     (default both)
#   TOL            relative tolerance for -c compare mode (default 0, i.e. exact)
#   VERIFY_MODE    light|full|none (default light). "light" runs the model
#                  NRUNS times in hash mode and compares compact full-field
#                  hashes, avoiding large baseline files. "full" preserves
#                  the original generate-baseline plus compare workflow.
#   NRUNS          number of correctness runs to perform, to check
#                  run-to-run reproducibility (default 2)
#   REPEAT         if >0, also do a dedicated timing run with this many
#                  repetitions (P3F "hot" iterations averaged); this uses the
#                  test's built-in -r timing mode, which reports the P3
#                  kernel's own internal timer (std::chrono::steady_clock
#                  wrapped tightly around the Kokkos::parallel_for compute
#                  kernel in p3_main_impl.hpp, excluding init) (default 20)
#   RESULTS_DIR    where to store logs/baselines/summary
#                  (default: a machine-specific <repo-root>/*-p3_test_results)
#   PROGRESS_INTERVAL  seconds between live heartbeat lines while a long
#                  timing executable invocation is still running (default 60)
#   NPAR           number of copies of this exact same test to launch
#                  concurrently on the compute node (default 1, i.e. no
#                  parallelism). Meant to exercise the node's parallel
#                  resources (4 GPUs on a GPU node, ~128 cores on a
#                  pm-cpu node) by running many identical instances at once
#                  and confirming they all produce the same fingerprint --
#                  a cheap way to check for e.g. race conditions or
#                  node-to-node variability, not just run-to-run (serial)
#                  reproducibility. Each instance is pinned to its own GPU
#                  (CUDA_VISIBLE_DEVICES, round-robin over GPUs) for GPU
#                  tags, or its own CPU core range (taskset, round-robin over
#                  CPUS_PER_INSTANCE cores) otherwise.
#   CPUS_PER_INSTANCE  cores given to each parallel instance in CPU mode
#                  (default: 128/NPAR, i.e. spread evenly across a pm-cpu
#                  node's 128 cores)
#   PROFILE        off|nsys|ncu (default off). GPU tags only; ignored on CPU.
#                  Wraps only the dedicated timing run (REPEAT reps) with a
#                  profiler and writes the report
#                  under RUN_DIR/profile.*:
#                    nsys - Nsight Systems timeline/API/kernel-duration
#                           trace (profile.nsys-rep + profile.txt stats).
#                           Good for overall time breakdown (H2D/D2H,
#                           kernel launch overhead, memory ops), but does
#                           NOT report FP32 vs FP64 instruction counts.
#                    ncu  - Nsight Compute, same P3 kernels/metrics as
#                           ./measure.sh, including per-kernel FP32
#                           (fadd/fmul/ffma) and FP64 (dadd/dmul/dfma)
#                           instruction counts (profile.ncu-rep +
#                           profile.txt). This is the tool that actually
#                           answers "how many single vs double precision
#                           ops did the kernel execute" -- nsys can't.
#                           Adds noticeable overhead (instrumented replay),
#                           so use a smaller REPEAT/NCOL than a normal
#                           timing run. With NPAR>1, only inst0 actually
#                           runs the ncu profiler (the other NPAR-1
#                           instances still do the timing/correctness
#                           parts, just skip profiling) since the
#                           per-kernel FP32/FP64 counts are deterministic
#                           across identical GPU instances and profiling
#                           all of them just adds redundant wall time (see
#                           notes.ncu.txt). Use NCU_LAUNCH_COUNT (and
#                           optionally NCU_LAUNCH_SKIP) env vars, passed
#                           through to ./measure.sh, to sample only the
#                           first N kernel launches instead of every launch
#                           across all NSTEPS -- a full unrestricted run
#                           needs ~3000 launches and does not fit in a
#                           30-min debug-QOS job.
#   NCU_LAUNCH_SKIP        CUDA kernel launches to skip when PROFILE=ncu
#                          (default 0).
#   NCU_LAUNCH_COUNT       CUDA kernel launches to collect when PROFILE=ncu
#                          (default 0, meaning omit the option/all launches).
#   NCU_KERNEL_NAME        ncu kernel-name filter (default: all P3 dispatch kernels).
#
# Live progress: timestamped "[HH:MM:SS] [<instance>] <message>" lines are
# printed to stdout as each stage starts/finishes (timing run, each
# correctness iteration, fingerprint, done), including elapsed seconds for
# each stage. Since sbatch scripts redirect stdout to the -o file (e.g.
# <machine>-p3_test_results/o.p3_test_<tag>.<jobid>.txt), this lets you `tail -f` that
# file to monitor an in-flight run instead of seeing nothing until the job
# finishes. With NPAR>1, every concurrent instance prints its own
# progress lines (labeled by run-dir basename, e.g. "inst0"), interleaved.
#
# Timing data recorded per run (under <tag>/<timestamp>[/inst<N>]/):
#   timing.log      P3 kernel-only time (internal steady_clock timer),
#                    averaged over REPEAT reps, excluding the cold first call
#   wallclock.log    wall-clock start/end/elapsed for every single exe
#                    invocation (timing run + each of the NRUNS correctness
#                    runs), covering full process launch + I/O + kernel, not
#                    just the kernel
#   summary.txt      both of the above plus knobs, fingerprint, pass/fail
# When NPAR>1, a top-level parallel_summary.txt additionally records whether
# every instance's fingerprint matched.
#
# Bigger/harder problem: increase NCOL, NLEV, NSTEPS, or REPEAT.
# Smaller/easier/faster problem: decrease them.
#
set -e

TEST_EXE="$1"
TAG="${2:-p3_test}"

if [ -z "$TEST_EXE" ] || [ ! -x "$TEST_EXE" ]; then
  echo "Error: run_p3_test.sh requires a valid path to the p3_run_and_cmp executable as \$1" >&2
  exit 1
fi
TEST_EXE="$(readlink -f "$TEST_EXE")"

NCOL=${NCOL:-3}
NLEV=${NLEV:-72}
NSTEPS=${NSTEPS:-6}
WARMUP=${WARMUP:-3}
DT=${DT:-300}
PREDICT_NC=${PREDICT_NC:-both}
PRESCRIBED_CCN=${PRESCRIBED_CCN:-both}
TOL=${TOL:-0}
VERIFY_MODE=${VERIFY_MODE:-light}
NRUNS=${NRUNS:-2}
REPEAT=${REPEAT:-20}
NPAR=${NPAR:-1}
PROFILE=${PROFILE:-off}
NCU_LAUNCH_SKIP=${NCU_LAUNCH_SKIP:-0}
NCU_LAUNCH_COUNT=${NCU_LAUNCH_COUNT:-0}
NCU_KERNEL_NAME=${NCU_KERNEL_NAME:-'regex:.*(p3_check_values|p3_main_init|p3_main_part1|p3_main_part2|p3_main_part3|cloud_sedimentation|rain_sedimentation|ice_sedimentation|homogeneous_freezing)_disp.*'}
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
case "$PROFILE" in
  off|nsys|ncu) ;;
  *)
    echo "Error: PROFILE must be off, nsys, or ncu; got '${PROFILE}'" >&2
    exit 1
    ;;
esac
case "$TAG" in
  *alvarez-gpu*) RESULTS_ROOT="${REPO_ROOT}/alvarez-gpu-p3_test_results" ;;
  *alvarez-cpu*) RESULTS_ROOT="${REPO_ROOT}/alvarez-cpu-p3_test_results" ;;
  *pm-cpu*) RESULTS_ROOT="${REPO_ROOT}/pm-cpu-p3_test_results" ;;
  *) RESULTS_ROOT="${REPO_ROOT}/pm-gpu-p3_test_results" ;;
esac
RESULTS_DIR=${RESULTS_DIR:-"${RESULTS_ROOT}/${TAG}"}

is_gpu_tag () {
  case "$1" in
    *gpu*|vista-gh*) return 0 ;;
    *) return 1 ;;
  esac
}

STAMP=$(date +%Y%m%d_%H%M%S)
TOP_RUN_DIR="${RESULTS_DIR}/${STAMP}"
mkdir -p "$TOP_RUN_DIR"

COMMON_ARGS=(-s "$NSTEPS" -w "$WARMUP" -dt "$DT" -i "$NCOL" -k "$NLEV"
             --predict-nc "$PREDICT_NC" --prescribed-ccn "$PRESCRIBED_CCN")

# Runs one full instance of the test (timing run + NRUNS correctness runs +
# fingerprint) into $1 (a run directory), optionally prefixing every exe
# invocation with a pinning command in $2 (e.g. "taskset -c 0-31" or
# "env CUDA_VISIBLE_DEVICES=2"). Writes RUN_DIR/summary.txt and, if a
# baseline file is produced, RUN_DIR/fingerprint.txt (just the hash, for easy
# cross-instance comparison).
run_one_instance () {
  local RUN_DIR="$1"; shift
  local PIN=("$@")
  local BASELINE_DIR="${RUN_DIR}/baselines"
  mkdir -p "$BASELINE_DIR"

  # Instance label used to prefix live-progress messages printed to stdout
  # (which lands in the sbatch -o file, e.g. o.p3_test_*.<jobid>.txt), so a
  # user tailing/reading that file while the job is still running can see
  # what's happening instead of a blank file until the very end.
  local INST_LABEL
  INST_LABEL="$(basename "$RUN_DIR")"
  progress () {
    printf '[%s] [%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$INST_LABEL" "$*"
  }

  local WALLCLOCK_LOG="${RUN_DIR}/wallclock.log"
  : > "$WALLCLOCK_LOG"
  # Runs the exe (with pinning prefix) with args, timing the whole process
  # (launch + I/O + kernel, not just the P3 kernel itself), and records
  # start/end/elapsed to wallclock.log. This is independent of/in addition
  # to the test's own internal std::chrono::steady_clock kernel timer (used
  # by -r/REPEAT above).
  run_timed () {
    local label="$1"; shift
    local t0 t1
    t0=$(date +%s.%N)
    "${PIN[@]}" "$@"
    local rc=$?
    t1=$(date +%s.%N)
    printf '%-20s start=%s end=%s elapsed=%.3fs\n' \
      "$label" "$t0" "$t1" "$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')" >> "$WALLCLOCK_LOG"
    return $rc
  }

  # Same as run_timed, but captures stdout/stderr to a log file and prints
  # periodic status lines while the command is still running. The P3
  # executable only emits the kernel timing summary after the -r timing run
  # finishes, so this keeps the sbatch output from looking stuck.
  run_timed_to_log_with_heartbeat () {
    local label="$1"; shift
    local log_file="$1"; shift
    local heartbeat_msg="$1"; shift
    local interval="${PROGRESS_INTERVAL:-60}"
    local t0 t1 hb_pid rc
    t0=$(date +%s.%N)
    (
      while :; do
        sleep "$interval" || exit 0
        local now elapsed
        now=$(date +%s.%N)
        elapsed=$(awk -v a="$t0" -v b="$now" 'BEGIN{printf "%.1f", b-a}')
        progress "${heartbeat_msg}; elapsed=${elapsed}s; log=${log_file}"
      done
    ) &
    hb_pid=$!
    set +e
    "${PIN[@]}" "$@" > "$log_file" 2>&1
    rc=$?
    set -e
    kill "$hb_pid" >/dev/null 2>&1 || true
    wait "$hb_pid" >/dev/null 2>&1 || true
    t1=$(date +%s.%N)
    printf '%-20s start=%s end=%s elapsed=%.3fs\n' \
      "$label" "$t0" "$t1" "$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')" >> "$WALLCLOCK_LOG"
    return $rc
  }

  local SUMMARY="${RUN_DIR}/summary.txt"
  {
    echo "P3 run_and_cmp test summary"
    echo "  tag            = ${TAG}"
    echo "  exe            = ${TEST_EXE}"
    echo "  run_dir        = ${RUN_DIR}"
    echo "  pinning        = ${PIN[*]:-(none)}"
    echo "  ncol           = ${NCOL}"
    echo "  nlev           = ${NLEV}"
    echo "  nsteps         = ${NSTEPS}"
    echo "  warmup         = ${WARMUP}"
    echo "  dt             = ${DT}"
    echo "  predict_nc     = ${PREDICT_NC}"
    echo "  prescribed_ccn = ${PRESCRIBED_CCN}"
    echo "  tol            = ${TOL}"
    echo "  verify_mode    = ${VERIFY_MODE}"
    echo "  nruns          = ${NRUNS}"
    echo "  repeat (perf)  = ${REPEAT}"
    echo
  } > "$SUMMARY"

  progress "starting: ncol=${NCOL} nlev=${NLEV} warmup=${WARMUP} nsteps=${NSTEPS} dt=${DT} nruns=${NRUNS} repeat=${REPEAT} pin=${PIN[*]:-(none)}"

  # --- 1. Timing run (optional, off if REPEAT<=0) ---
  if [ "$REPEAT" -gt 0 ]; then
    echo "==========================================================" >> "$SUMMARY"
    echo " Timing run (repeat=${REPEAT}, no I/O)" >> "$SUMMARY"
    echo "==========================================================" >> "$SUMMARY"
    local PERF_LOG="${RUN_DIR}/timing.log"
    progress "timing run starting (repeat=${REPEAT})..."
    local TIMING_T0 TIMING_T1
    TIMING_T0=$(date +%s.%N)
    run_timed_to_log_with_heartbeat "timing_run" "$PERF_LOG" \
      "timing run still running (repeat=${REPEAT}, ncol=${NCOL}, nlev=${NLEV}, nsteps=${NSTEPS})" \
      "$TEST_EXE" -g -b "$BASELINE_DIR" "${COMMON_ARGS[@]}" -r "$REPEAT"
    TIMING_T1=$(date +%s.%N)
    echo "--- Timing results (Time = seconds/call, averaged over ${REPEAT} reps, per case) ---" >> "$SUMMARY"
    grep -E "^Running P3|^Timing summary:|^Time = |^(Measured steps|First measured step|Step [0-9]+):" "$PERF_LOG" >> "$SUMMARY"
    echo >> "$SUMMARY"
    progress "timing run done in $(awk -v a="$TIMING_T0" -v b="$TIMING_T1" 'BEGIN{printf "%.1f", b-a}')s; kernel timing follows"
    while IFS= read -r line; do
      progress "timing: ${line}"
    done < <(grep -E "^Running P3|^Timing summary:|^Time = |^(Measured steps|First measured step|Step [0-9]+):" "$PERF_LOG")
  fi

  # --- 1b. Optional profiling run (nsys or ncu), GPU only ---
  if [ "$PROFILE" != "off" ]; then
    if is_gpu_tag "$TAG"; then
        local PROFILE_TXT="${RUN_DIR}/profile.txt"
        if [ "$PROFILE" = "nsys" ]; then
          progress "nsys profile starting..."
          "${PIN[@]}" nsys profile -o "${RUN_DIR}/profile" --force-overwrite=true \
            --stats=true "$TEST_EXE" -g -b "$BASELINE_DIR" "${COMMON_ARGS[@]}" -r 5 \
            > "$PROFILE_TXT" 2>&1 || true
          progress "nsys profile done -> ${RUN_DIR}/profile.nsys-rep (+ profile.txt)"
        elif [ "$PROFILE" = "ncu" ]; then
          # ncu per-kernel FP32/FP64 instruction counts are deterministic
          # (same kernel code, same problem size), so profiling all NPAR
          # identical GPU instances gives zero new information -- just NPARx
          # the redundant profile.txt files and GPU time, which is exactly
          # what makes runs blow past the debug-QOS 30-min walltime (see
          # notes.ncu.txt). NPAR is still kept for the correctness/timing
          # race-condition check above; only inst0 (or the sole instance
          # when NPAR<=1) actually runs the ncu profiler.
          if [ "$NPAR" -gt 1 ] && [ "$INST_LABEL" != "inst0" ]; then
            progress "skipping ncu profile on ${INST_LABEL} (NPAR>1; only inst0 profiles, see notes.ncu.txt)"
          else
            progress "ncu profile starting (FP32/FP64 instruction counts)..."
            local NCU_ARGS=(--set full --target-processes all
                            --kernel-name-base demangled
                            --launch-skip "$NCU_LAUNCH_SKIP"
                            -o "${RUN_DIR}/profile")
            if [ "$NCU_LAUNCH_COUNT" -gt 0 ]; then
              NCU_ARGS+=(--launch-count "$NCU_LAUNCH_COUNT")
            fi
            if [ -n "$NCU_KERNEL_NAME" ]; then
              NCU_ARGS+=(--kernel-name "$NCU_KERNEL_NAME")
            fi
            "${PIN[@]}" ncu "${NCU_ARGS[@]}" "$TEST_EXE" -g -b "$BASELINE_DIR" \
              "${COMMON_ARGS[@]}" -r 5 > "$PROFILE_TXT" 2>&1
            progress "ncu profile done -> ${RUN_DIR}/profile.ncu-rep (+ profile.txt)"
          fi
        else
          progress "unknown PROFILE=${PROFILE}, skipping (expected nsys|ncu|off)"
        fi
    else
        progress "PROFILE=${PROFILE} requested but tag=${TAG} is not GPU; skipping"
    fi
  fi

  # --- 2. Correctness / reproducibility runs ---
  echo "==========================================================" >> "$SUMMARY"
  echo " Correctness/reproducibility runs (mode=${VERIFY_MODE}, n=${NRUNS}, tol=${TOL})" >> "$SUMMARY"
  echo "==========================================================" >> "$SUMMARY"

  local i LOG HASH FIRST_HASH
  if [ "$VERIFY_MODE" = "light" ]; then
    for i in $(seq 1 "$NRUNS"); do
      LOG="${RUN_DIR}/hash_iter${i}.log"
      local ITER_T0 ITER_T1
      ITER_T0=$(date +%s.%N)
      progress "hash check run ${i}/${NRUNS} starting..."
      run_timed "hash_iter${i}" "$TEST_EXE" --hash "${COMMON_ARGS[@]}" > "$LOG" 2>&1
      ITER_T1=$(date +%s.%N)
      HASH=$(grep '^P3 hash ' "$LOG" | sha256sum | awk '{print $1}')
      echo "hash_iter${i} ${HASH}" >> "${RUN_DIR}/hashes.txt"
      if [ "$i" -eq 1 ]; then
        FIRST_HASH="$HASH"
        echo "$HASH" > "${RUN_DIR}/fingerprint.txt"
      elif [ "$HASH" != "$FIRST_HASH" ]; then
        progress "hash check run ${i}/${NRUNS} mismatch: ${HASH} != ${FIRST_HASH}"
        return 1
      fi
      progress "hash check run ${i}/${NRUNS} done in $(awk -v a="$ITER_T0" -v b="$ITER_T1" 'BEGIN{printf "%.1f", b-a}')s; fingerprint=${HASH}"
    done
    echo "Light hash fingerprints:" >> "$SUMMARY"
    cat "${RUN_DIR}/hashes.txt" >> "$SUMMARY"
  elif [ "$VERIFY_MODE" = "full" ]; then
    # Run 1 generates a reference baseline (full field dump). Runs 2..NRUNS
    # compare against that same baseline at $TOL.
    for i in $(seq 1 "$NRUNS"); do
      LOG="${RUN_DIR}/run_iter${i}.log"
      local ITER_T0 ITER_T1
      ITER_T0=$(date +%s.%N)
      if [ "$i" -eq 1 ]; then
        progress "correctness run ${i}/${NRUNS} starting (generate baseline)..."
        run_timed "run_iter1_gen" "$TEST_EXE" -g -b "$BASELINE_DIR" "${COMMON_ARGS[@]}" > "$LOG" 2>&1
      else
        progress "correctness run ${i}/${NRUNS} starting (compare, tol=${TOL})..."
        run_timed "run_iter${i}_cmp" "$TEST_EXE" -c -b "$BASELINE_DIR" -t "$TOL" "${COMMON_ARGS[@]}" > "$LOG" 2>&1
      fi
      ITER_T1=$(date +%s.%N)
      progress "correctness run ${i}/${NRUNS} done in $(awk -v a="$ITER_T0" -v b="$ITER_T1" 'BEGIN{printf "%.1f", b-a}')s"
    done

    local BASELINE_FILE
    BASELINE_FILE=$(find "$BASELINE_DIR" -maxdepth 1 -name 'p3_run_and_cmp.baseline*' | head -1)
    if [ -n "$BASELINE_FILE" ]; then
      HASH=$(sha256sum "$BASELINE_FILE" | awk '{print $1}')
      echo >> "$SUMMARY"
      echo "Result fingerprint (sha256 of full field output, ${BASELINE_FILE##*/}):" >> "$SUMMARY"
      echo "  ${HASH}" >> "$SUMMARY"
      echo "$HASH" > "${RUN_DIR}/fingerprint.txt"
      progress "fingerprint = ${HASH}"
    fi
  elif [ "$VERIFY_MODE" = "none" ]; then
    progress "VERIFY_MODE=none; skipping correctness checks"
  else
    progress "unknown VERIFY_MODE=${VERIFY_MODE} (expected light|full|none)"
    return 1
  fi

  echo >> "$SUMMARY"
  echo "--- Wall-clock time per exe invocation (process launch + I/O + kernel) ---" >> "$SUMMARY"
  cat "$WALLCLOCK_LOG" >> "$SUMMARY"

  echo >> "$SUMMARY"
  echo "--- Per-run pass/fail (nerr from compare, 0 == identical within tol) ---" >> "$SUMMARY"
  grep -H -E "Ref impl failed|^Comparing with|^Running with|^P3 hash " "${RUN_DIR}"/run_iter*.log "${RUN_DIR}"/hash_iter*.log >> "$SUMMARY" 2>/dev/null || true

  echo >> "$SUMMARY"
  echo "Results saved under: ${RUN_DIR}" >> "$SUMMARY"
  progress "finished"
}

# --- DCGM profiling counters conflict with ncu's exclusive access to the
# same hardware performance counters. On Perlmutter (and similar systems),
# DCGM's background health-monitoring daemon runs on every GPU node and
# holds those counters, so `ncu` fails with "Profiling failed because a
# driver resource was unavailable" / "Failed to profile
# query_cuda_kernel_arch" unless DCGM is paused for the duration of the ncu
# run. Pause it once here (covers both the single-instance and the NPAR>1
# concurrent-instances cases below, since all instances on this node share
# the same DCGM daemon) and resume it unconditionally on exit.
DCGM_PAUSED=0
DCGM_REPAUSE_PID=""
if is_gpu_tag "$TAG"; then
    if [ "$PROFILE" = "ncu" ]; then
      if command -v dcgmi >/dev/null 2>&1; then
        echo "Pausing DCGM profiling counters (required for ncu to get exclusive access to perf counters)..."
        if dcgmi profile --pause >/dev/null 2>&1; then
          DCGM_PAUSED=1
        else
          echo "WARNING: 'dcgmi profile --pause' failed; ncu profiling may fail with 'driver resource unavailable'." >&2
        fi
      else
        echo "WARNING: dcgmi not found; ncu profiling may fail with 'driver resource unavailable' if DCGM is active on this node." >&2
      fi

      # Something on Perlmutter (a periodic node-health probe/cron) appears
      # to re-enable DCGM mid-job independent of our one-shot pause above
      # (observed: all concurrent ncu profiling passes across multiple nodes
      # froze simultaneously at a suspiciously round wall-clock minute
      # boundary, see notes.ncu.txt). Since a resumed DCGM makes `ncu` hang
      # (rather than fail fast like the initial "driver resource
      # unavailable" error), guard against it by re-issuing
      # `dcgmi profile --pause` in the background every few minutes for the
      # duration of this script. This is a no-op if DCGM is already paused.
      if [ "$DCGM_PAUSED" -eq 1 ]; then
        (
          while sleep 120; do
            dcgmi profile --pause >/dev/null 2>&1 || true
          done
        ) &
        DCGM_REPAUSE_PID=$!
      fi
    fi
fi
resume_dcgm() {
  if [ -n "$DCGM_REPAUSE_PID" ]; then
    kill "$DCGM_REPAUSE_PID" >/dev/null 2>&1 || true
    wait "$DCGM_REPAUSE_PID" 2>/dev/null || true
    DCGM_REPAUSE_PID=""
  fi
  if [ "$DCGM_PAUSED" -eq 1 ]; then
    dcgmi profile --resume >/dev/null 2>&1 || echo "WARNING: 'dcgmi profile --resume' failed" >&2
    DCGM_PAUSED=0
  fi
}
trap resume_dcgm EXIT

if [ "$NPAR" -le 1 ]; then
  # --- Single instance: original behavior, streaming to stdout as we go ---
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting single-instance run (tag=${TAG}) under ${TOP_RUN_DIR}"
  run_one_instance "$TOP_RUN_DIR"
  cat "${TOP_RUN_DIR}/summary.txt"
else
  # --- Multiple concurrent instances of the exact same run ---
  # This is meant to make use of a compute node's parallel resources (4 GPUs
  # on a GPU node, ~128 cores on a CPU node) by firing off NPAR identical copies at
  # once and confirming every one produces the same fingerprint, i.e. a
  # cheap concurrency/node-variability check on top of the serial NRUNS
  # reproducibility check above.
  IS_GPU=0
  if is_gpu_tag "$TAG"; then IS_GPU=1; fi

  CPUS_PER_INSTANCE=${CPUS_PER_INSTANCE:-$((128 / NPAR))}
  if [ "$CPUS_PER_INSTANCE" -lt 1 ]; then CPUS_PER_INSTANCE=1; fi

  echo "Launching ${NPAR} concurrent instances (tag=${TAG}, gpu_mode=${IS_GPU}) under ${TOP_RUN_DIR}"
  PIDS=()
  for n in $(seq 0 $((NPAR - 1))); do
    INST_DIR="${TOP_RUN_DIR}/inst${n}"
    mkdir -p "$INST_DIR"
    if [ "$IS_GPU" -eq 1 ]; then
      GPU_COUNT=${GPU_COUNT:-4}
      GPU_ID=$((n % GPU_COUNT))
      PIN=(env CUDA_VISIBLE_DEVICES="$GPU_ID" HIP_VISIBLE_DEVICES="$GPU_ID")
    else
      LO=$((n * CPUS_PER_INSTANCE))
      HI=$((LO + CPUS_PER_INSTANCE - 1))
      PIN=(taskset -c "${LO}-${HI}")
    fi
    ( run_one_instance "$INST_DIR" "${PIN[@]}" ) &
    PIDS+=($!)
  done

  FAIL=0
  for pid in "${PIDS[@]}"; do
    wait "$pid" || FAIL=1
  done

  PAR_SUMMARY="${TOP_RUN_DIR}/parallel_summary.txt"
  {
    echo "Parallel P3 run_and_cmp summary (tag=${TAG}, NPAR=${NPAR})"
    echo
    for n in $(seq 0 $((NPAR - 1))); do
      FP_FILE="${TOP_RUN_DIR}/inst${n}/fingerprint.txt"
      if [ -f "$FP_FILE" ]; then
        printf "  inst%-3d fingerprint = %s\n" "$n" "$(cat "$FP_FILE")"
      else
        printf "  inst%-3d fingerprint = (missing)\n" "$n"
      fi
    done
    echo
    if [ "$VERIFY_MODE" = "none" ]; then
      echo "PASS: VERIFY_MODE=none; no fingerprints were expected."
    else
      NFP=$(cat "${TOP_RUN_DIR}"/inst*/fingerprint.txt 2>/dev/null | wc -l)
      UNIQ=$(cat "${TOP_RUN_DIR}"/inst*/fingerprint.txt 2>/dev/null | sort -u | wc -l)
      if [ "$UNIQ" -eq 1 ] && [ "$NFP" -eq "$NPAR" ]; then
        echo "PASS: all ${NPAR} concurrent instances produced identical fingerprints."
      else
        echo "FAIL: got ${NFP}/${NPAR} fingerprints with ${UNIQ} distinct value(s) (expected ${NPAR}/${NPAR}, 1 distinct)."
      fi
    fi
  } | tee "$PAR_SUMMARY"

  if [ "$FAIL" -ne 0 ]; then
    echo "Warning: one or more instances exited non-zero; check inst*/summary.txt" >&2
  fi

  echo "Results saved under: ${TOP_RUN_DIR}"
fi
