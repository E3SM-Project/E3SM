#!/bin/bash
#SBATCH -A e3sm
#SBATCH -q debug
#SBATCH -C gpu
## #SBATCH -C  gpu&hbm40g
#SBATCH -N 1
#SBATCH --gpus=4
#SBATCH -t 00:30:00
#SBATCH -J p3_test_gpu-sp
#SBATCH -o alvarez-gpu-p3_test_results/o.p3_test_alvarez-gpu-sp.%j.txt
#
# Ready-to-submit batch job for the Alvarez GPU, SP p3_run_and_cmp test.
# Requests 1 GPU node (4 GPUs) for 10 minutes on the debug QOS.
#
# Build the executable on the login node FIRST (not part of this job):
#   ./build_p3_test_alvarez-gpu-sp.sh
#
# Then submit this script from the repo root:
#   sbatch sbatch_p3_test_alvarez-gpu-sp.sh
#
# Edit the knobs below and just `sbatch sbatch_p3_test_alvarez-gpu-sp.sh` --
# no need to set env vars on the command line.
# See scripts_p3_test/run_p3_test.sh for the full knob list.
set -e
export MACHINE_NAME="${MACHINE_NAME:-alvarez-gpu}"
export NPAR=1      # number of concurrent instances (one GPU per instance)
export NCOL=20000      # number of columns (problem size knob)
export NLEV=128     # number of vertical levels (problem size knob)
export NSTEPS=6     # number of timesteps to run
export DT=300       # timestep length in seconds
export WARMUP=3     # warmup timesteps excluded from measured timing
export REPEAT=20    # number of "hot" repetitions in the dedicated timing run (kernel-only timer)
export VERIFY_MODE=none # light=hash check, full=baseline file compare, none=timing only
export NRUNS=2      # number of correctness runs to check run-to-run reproducibility
export TOL=0        # relative tolerance for correctness comparison (0 = exact/bit-for-bit)
# Profiling (GPU only, off by default -- comment out to run without it):
#   nsys = Nsight Systems timeline/kernel-duration trace (profile.nsys-rep)
#   ncu  = Nsight Compute, incl. FP32 vs FP64 instruction counts (profile.ncu-rep)
#          (this is the one that answers "how much single vs double precision
#          math", nsys does not report instruction mix)
#export PROFILE=ncu
# Skip the ~118 Kokkos-startup/View-construction launches first (confirmed
# against out-ncu-dp.txt: all 8 distinct physics disp kernels first appear
# at launch indices 119-126), then profile a 60-launch window -- comfortably
# covering p3_main_init/part1/part2/part3_disp, {cloud,rain,ice}_sedimentation_disp,
# and homogeneous_freezing_disp -- instead of profiling all ~3000 launches
# across every timestep/rep (FP32/FP64 instruction mix per kernel type is
# deterministic, so one representative window is enough). This cuts wall
# time to a small fraction of the ~42-45 min a full profile needs,
# comfortably inside the 30-min debug QOS cap. See notes.ncu.txt.
export NCU_LAUNCH_SKIP=100
export NCU_LAUNCH_COUNT=60

# Alvarez does not provide the Environment Modules command.  The executable
# was built with the CUDA 12.9 toolkit, so add its runtime library explicitly.
CUDA_LIB_DIR="${CUDA_LIB_DIR:-/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda/12.9/lib64}"
export LD_LIBRARY_PATH="${CUDA_LIB_DIR}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

echo "hostname=$(hostname)"
echo "CUDA_LIB_DIR=${CUDA_LIB_DIR}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
ldd "${SLURM_SUBMIT_DIR}/build_p3_test_alvarez-gpu-SP/src/physics/p3/tests/p3_run_and_cmp" | grep -E 'cuda|not found' || true

# NOTE: Slurm copies this script into a spool directory before running it,
# so BASH_SOURCE[0]'s dirname is NOT the repo root at run time. Use
# SLURM_SUBMIT_DIR (the directory `sbatch` was invoked from) instead.
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
"${REPO_ROOT}/run_p3_test_alvarez-gpu-sp.sh" "$MACHINE_NAME"
