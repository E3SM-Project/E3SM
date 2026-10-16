#!/bin/bash

set -e

usage() {
  echo "Usage: $0 [--compiler <name>] [--driver <name>] [--din-loc-root <path>] [--case-tag <tag>] [--nsegspc <val>] [--hist-nhtfrq <val>] [--hist-mfilt <val>]"
  exit 1
}

RES=t-lr25_t-lr25
COMPSET=I1850GSWCNPTOP
COMPILER=intel
DRIVER=moab
#DIN_LOC_ROOT="${SCRATCH}/inputdata"
DIN_LOC_ROOT="/dvs_ro/cfs/cdirs/e3sm/inputdata"
CASE_TAG="n008.exp5_cftprune"
NSEGSPC=""
HIST_NHTFRQ=""
HIST_MFILT="1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --compiler)
      [[ $# -ge 2 ]] || usage
      COMPILER=$2
      shift 2
      ;;
    --driver)
      [[ $# -ge 2 ]] || usage
      DRIVER=$2
      shift 2
      ;;
    --din-loc-root)
      [[ $# -ge 2 ]] || usage
      DIN_LOC_ROOT=$2
      shift 2
      ;;
    --case-tag)
      [[ $# -ge 2 ]] || usage
      CASE_TAG=$2
      shift 2
      ;;
    --nsegspc)
      [[ $# -ge 2 ]] || usage
      NSEGSPC=$2
      shift 2
      ;;
    --hist-nhtfrq)
      [[ $# -ge 2 ]] || usage
      HIST_NHTFRQ=$2
      shift 2
      ;;
    --hist-mfilt)
      [[ $# -ge 2 ]] || usage
      HIST_MFILT=$2
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "ERROR: Unknown argument '$1'"
      usage
      ;;
  esac
done

NTASKS=1024

MACH=pm-cpu
CHECKOUT=bishtgautam_lnd_add-t-meshes

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC_DIR=${SCRIPT_DIR}
SCRIPT_PATH=${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")

cd ${SRC_DIR}
GIT_HASH=`git log -n 1 --format=%h`
CASE_NAME="nocrop.${RES}.${COMPSET}.${COMPILER}.${CASE_TAG}"
CASE_ROOT="${SCRATCH}/e3sm_scratch/${MACH}/${CHECKOUT}/${CASE_NAME}"
CASE_BUILD_DIR=${CASE_ROOT}/build
CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
CASE_RUN_DIR=${CASE_ROOT}/run

cd ${SRC_DIR}/cime/scripts

COMPILER_OPT=()
[[ -n ${COMPILER} ]] && COMPILER_OPT=(--compiler "${COMPILER}")
DRIVER_OPT=()
[[ -n ${DRIVER} ]] && DRIVER_OPT=(--driver "${DRIVER}")

./create_newcase --case ${CASE_NAME} \
--output-root ${CASE_ROOT} \
--script-root ${CASE_SCRIPTS_DIR} \
--handle-preexisting-dirs u \
--res ${RES} --machine ${MACH} "${COMPILER_OPT[@]}" "${DRIVER_OPT[@]}" --compset ${COMPSET}

SCRIPT_PROVENANCE_DIR=${CASE_SCRIPTS_DIR}/run_script_provenance
mkdir -p "${SCRIPT_PROVENANCE_DIR}"
SCRIPT_PROVENANCE_NAME=$(basename "${SCRIPT_PATH}").$(date +%Y%m%d-%H%M%S)
cp -p "${SCRIPT_PATH}" "${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}"

cd ${CASE_SCRIPTS_DIR}

./xmlchange EXEROOT=${CASE_BUILD_DIR}
./xmlchange RUNDIR=${CASE_RUN_DIR}
./xmlchange DIN_LOC_ROOT="${DIN_LOC_ROOT}"
./xmlchange DIN_LOC_ROOT_CLMFORC="${DIN_LOC_ROOT}/atm/datm7"

./xmlchange NTASKS=${NTASKS}
./xmlchange REST_N=10,REST_OPTION=never

./case.setup

if [[ -n "${HIST_NHTFRQ}" ]]; then
  echo "hist_nhtfrq = ${HIST_NHTFRQ}" >> user_nl_elm
  echo "hist_mfilt = ${HIST_MFILT}" >> user_nl_elm
fi
if [[ -n "${NSEGSPC}" ]]; then
  echo "nsegspc = ${NSEGSPC}" >> user_nl_elm
fi

./case.build


chmod -R og+rX ${CASE_ROOT}

./case.submit -a="--qos=regular"
