#!/bin/bash
# Submit a p3_run_and_cmp test run as a batch job on Perlmutter (NERSC) or
# Vista (TACC).
#
# This only submits the *run* step (see run_p3_test_*.sh) -- build the
# executable first on the login node with the matching build_p3_test_*.sh,
# since compiling on a compute-node allocation just wastes reserved
# node-time.
#
# Usage:
#   scripts_p3_test/submit_p3_test.sh <tag> [extra env=val ...]
#
# <tag> is one of: pm-gpu-dp pm-gpu-sp alvarez-gpu-dp alvarez-gpu-sp
# pm-cpu-dp pm-cpu-sp vista-gh-dp vista-gh-sp
#
# Examples:
#   scripts_p3_test/submit_p3_test.sh pm-gpu-dp
#   scripts_p3_test/submit_p3_test.sh pm-gpu-dp NPAR=4 NCOL=50
#   scripts_p3_test/submit_p3_test.sh pm-gpu-sp PROFILE=ncu NCOL=50 NCU_LAUNCH_COUNT=20
#   scripts_p3_test/submit_p3_test.sh pm-cpu-dp NPAR=16 CPUS_PER_INSTANCE=8
#   scripts_p3_test/submit_p3_test.sh vista-gh-dp NCOL=50
#
# Edit ACCOUNT/QOS/PARTITION below, or export SBATCH_ACCOUNT/SBATCH_QOS/
# SBATCH_PARTITION before calling, to match your allocation.
set -e

TAG="$1"; shift || true
if [ -z "$TAG" ]; then
  echo "Usage: submit_p3_test.sh <tag: pm-gpu-dp|pm-gpu-sp|alvarez-gpu-dp|alvarez-gpu-sp|pm-cpu-dp|pm-cpu-sp|vista-gh-dp|vista-gh-sp> [ENV=val ...]" >&2
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_TAG="$TAG"
MACHINE_NAME="${TAG%-dp}"
MACHINE_NAME="${MACHINE_NAME%-sp}"
case "$TAG" in
  alvarez-gpu-dp|alvarez-gpu-sp) RUN_TAG="$TAG" ;;
esac
RUN_SCRIPT="${REPO_ROOT}/run_p3_test_${RUN_TAG}.sh"
if [ ! -x "$RUN_SCRIPT" ]; then
  echo "Error: no run script for tag '${TAG}' (looked for ${RUN_SCRIPT})" >&2
  exit 1
fi

TIME="${SBATCH_TIME:-00:30:00}"

SITE_OPTS=()
case "$TAG" in
  pm-gpu*|alvarez-gpu*)
    ACCOUNT="${SBATCH_ACCOUNT:-e3sm}"
    QOS="${SBATCH_QOS:-debug}"
    SITE_OPTS=(-A "$ACCOUNT" -q "$QOS" -C gpu --gpus=4)
    ;;
  pm-cpu*)
    ACCOUNT="${SBATCH_ACCOUNT:-e3sm}"
    QOS="${SBATCH_QOS:-debug}"
    SITE_OPTS=(-A "$ACCOUNT" -q "$QOS" -C cpu)
    ;;
  vista-gh*)
    ACCOUNT="${SBATCH_ACCOUNT:-CDA24017}"
    PARTITION="${SBATCH_PARTITION:-gh-dev}"
    SITE_OPTS=(-A "$ACCOUNT" -p "$PARTITION")
    ;;
  *) echo "Error: unsupported machine tag '${TAG}'" >&2; exit 1 ;;
esac

# Remaining args are ENV=val pairs (e.g. NPAR=4 NCOL=50) forwarded to the run
# script inside the job.
ENV_ARGS=("$@")

JOBSCRIPT=$(mktemp)
{
  echo "#!/bin/bash"
  echo "#SBATCH -N 1"
  echo "#SBATCH -t ${TIME}"
  echo "#SBATCH -J p3_test_${TAG}"
  case "$TAG" in
    *alvarez-gpu*) RESULTS_ROOT="${REPO_ROOT}/alvarez-gpu-p3_test_results" ;;
    *) RESULTS_ROOT="${REPO_ROOT}/pm-gpu-p3_test_results" ;;
  esac
  echo "#SBATCH -o ${RESULTS_ROOT}/o.p3_test_${TAG}.%j.txt"
  for kv in "${ENV_ARGS[@]}"; do
    echo "export ${kv}"
  done
  echo "export MACHINE_NAME=${MACHINE_NAME}"
  echo "\"${RUN_SCRIPT}\""
} > "$JOBSCRIPT"

mkdir -p "${RESULTS_ROOT}"
sbatch "${SITE_OPTS[@]}" "$JOBSCRIPT"
