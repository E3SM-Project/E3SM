#!/bin/bash -fe

# E3SM run script for FME (Full Model Emulation) production.
#
# Generates ACE/Samudra training data via the fme_output testmod, which
# configures online horizontal remapping, vertical coarsening, derived
# diagnostics, and column-integrated fields for EAM, MPAS-O, and MPAS-SI.
#
# Adapted from the standard E3SM run_e3sm.template.sh. Key differences:
#   * --user-mods-dirs points at cime_config/testmods_dirs/allactive/fme_output
#   * user_nl() is empty: the testmod is the single source of truth for
#     EAM/MPAS-O/MPAS-SI namelists -- including the FME AM enable/config and
#     the disable of non-FME AMs (globalStats, regionalStatistics,
#     timeSeriesStats*) that would otherwise add storage overhead. Add
#     fincl2/fincl3 here if you want extra diagnostic tapes alongside the
#     FME tape.
#   * Production restart cadence is 5 yr. With the 2026-05-01 restart fix
#     (append-mode reopen + accumulator sidecar + frame tracking +
#     compute_on_startup=.false.; see AGENTS.md gotchas #11 and #29) the
#     model is BFB across restart, so the cadence choice is now about
#     wallclock budget and output-file size rather than restart-clobber
#     mitigation. 5-yr cadence with monthly file rotation keeps each
#     restart lookup small while still landing on a year-boundary, which
#     means leg N+1 always starts a fresh `*.YYYY-01.nc` (no append path
#     exercised in production).

main() {

# --- Configuration flags ----
#
# Site-specific knobs (MACHINE, PROJECT, PELAYOUT, WALLTIME, RUN_REF*, and
# the scratch root via $PSCRATCH) are env-overridable. Defaults below are
# the pm-cpu/SamudrACE production combo. To run elsewhere, export e.g.
#   MACHINE=chrysalis PSCRATCH=/lcrc/group/e3sm/$USER/scratch ./run_e3sm.fme.sh
# COMPILER is intentionally omitted from create_newcase so each machine's
# default compiler is used; pin it explicitly here if you need a non-default.

readonly MACHINE=${MACHINE:-pm-cpu}
readonly PROJECT="${PROJECT:-$(sacctmgr show user $USER format=DefaultAccount | tail -n1 | tr -d ' ')}"

readonly COMPSET="WCYCL1850"
readonly RESOLUTION="ne30pg2_r05_IcoswISC30E3r5"
# BEFORE RUNNING: set CASE_BASE to a descriptive identifier. CASE_TAG is
# an optional dotted suffix (env-overridable) for launching several
# variants side-by-side without name collisions, e.g.
#   CASE_TAG=intel ./run_e3sm.fme.sh   -> v3.LR.piControl.aigo.intel
#   CASE_TAG=test1 ./run_e3sm.fme.sh   -> v3.LR.piControl.aigo.test1
readonly CASE_BASE="${CASE_BASE:-v3.LR.piControl.aigo}"
readonly CASE_TAG="${CASE_TAG:-}"
readonly CASE_NAME="${CASE_BASE}${CASE_TAG:+.${CASE_TAG}}"
# readonly CASE_GROUP="samudrace_v3"

# Code and compilation
# BEFORE RUNNING: set CHECKOUT to a date string like 20260430
readonly CHECKOUT="latest"
readonly BRANCH="mahf708/fme/aigo"
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options.
#
# We branch from a previously spun-up v3 LR piControl rather than cold-starting
# (cold start would discard ~50 yr of ocean spinup). 'branch' = exact restart,
# bit-identical model evolution, just with FME diagnostics added on top.
# Use 'hybrid' instead if the source case had different physics/forcing.
readonly MODEL_START_TYPE="hybrid"   # 'initial', 'continue', 'branch', 'hybrid'

# BEFORE RUNNING (branch/hybrid only): fill in the three RUN_REF* values
# below to match the spun-up piControl restart you're branching from.
# START_DATE should match RUN_REFDATE for a 'branch' run. All three are
# env-overridable; the PLACEHOLDER defaults trip the guard further down.
readonly RUN_REFDIR="${RUN_REFDIR:-/pscratch/sd/m/mahf708/v3.LR.piControl/archive/rest/0401-01-01-00000}"
readonly RUN_REFCASE="${RUN_REFCASE:-v3.LR.piControl}"
readonly RUN_REFDATE="${RUN_REFDATE:-0401-01-01}"
readonly START_DATE="${START_DATE:-${RUN_REFDATE}}"

# GET_REFCASE=TRUE asks CIME to copy refcase restart files from RUN_REFDIR.
# Set to FALSE if you've pre-staged the files into CASE_RUN_DIR yourself.
readonly GET_REFCASE=TRUE

# Set paths
readonly CASE_ROOT="${PSCRATCH}/E3SMv3/${CASE_NAME}"
# CODE_ROOT defaults to a self-contained checkout under CASE_ROOT so each
# case is self-describing (source + build + run + archive co-located, no
# ambiguity about which sha was actually compiled). Override by exporting
# CODE_ROOT=/path/to/shared/checkout to point several cases at one tree
# (saves ~5 GB of git history + submodules per case, at the cost of
# losing per-case provenance).
readonly CODE_ROOT="${CODE_ROOT:-${CASE_ROOT}/code/${CHECKOUT}}"

readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# FME testmod directory (applied via --user-mods-dirs at create_newcase time)
readonly FME_TESTMOD="${CODE_ROOT}/cime_config/testmods_dirs/allactive/fme_output"

# Define type of run
#  short tests: 'XS_2x5_ndays', 'S_1x10_ndays', etc. (same scheme as upstream)
#  or 'production' for the SamudrACE 100-yr piControl tape
# Override via RUN_LAYOUT env var (e.g. RUN_LAYOUT='XS_1x2_ndays').
readonly run="${RUN_LAYOUT:-production}"

if [ "${run}" != "production" ]; then
  echo "setting up Short test simulations: ${run}"
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}

  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/tests/${run}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/tests/${run}/run
  readonly PELAYOUT=${PELAYOUT:-${layout}}
  readonly WALLTIME=${WALLTIME:-2:00:00}
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false

else

  # Production: 100-yr SamudrACE training data run.
  # 5-yr segments × 20 = 100 yr. Year-boundary restart means each new leg
  # opens a fresh monthly file (`*.YYYY-01.remapped.nc`) with no append
  # logic exercised. The append-mode + sidecar machinery is in place
  # (see AGENTS.md #29) and keeps mid-window restarts safe, but the
  # production cadence keeps things simple by avoiding that path.
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  # PELAYOUT and WALLTIME defaults are pm-cpu-tuned. On sites with
  # shorter max walltime (e.g. 24h on chrysalis), export WALLTIME=24:00:00.
  readonly PELAYOUT=${PELAYOUT:-L}
  readonly WALLTIME=${WALLTIME:-34:00:00}
  readonly STOP_OPTION="nyears"
  readonly STOP_N="5"
  readonly REST_OPTION="nyears"
  readonly REST_N="5"
  readonly RESUBMIT="19"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history
readonly HIST_OPTION="nyears"
readonly HIST_N="5"

readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=true
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

# --- Now, do the work ---

umask 022

# Refuse to run with PLACEHOLDER refcase values still in place.
if [ "${MODEL_START_TYPE,,}" == "branch" ] || [ "${MODEL_START_TYPE,,}" == "hybrid" ]; then
    if [[ "${RUN_REFDIR}" == *PLACEHOLDER* ]] || \
       [[ "${RUN_REFCASE}" == *PLACEHOLDER* ]] || \
       [[ "${CASE_NAME}" == *PLACEHOLDER* ]]; then
        echo "ERROR: PLACEHOLDER refcase values still present. Edit RUN_REFDIR,"
        echo "       RUN_REFCASE, RUN_REFDATE, and CASE_NAME at the top of the"
        echo "       script before submitting."
        exit 10
    fi
fi

fetch_code
create_newcase
custom_pelayout
case_setup
case_build
runtime_options
copy_script
case_submit

echo $'\n----- All done -----\n'

}

# =======================
# Custom user_nl settings
# =======================
#
# The fme_output testmod (applied via --user-mods-dirs) is the single
# source of truth for EAM/MPAS-O/MPAS-SI namelists. It configures:
#   * FME analysis members enabled with daily averaging
#     (fmeDepthCoarsening, fmeDerivedFields, fmeSeaiceDerivedFields)
#   * fmeVerticalReduce disabled (not in SamudrACE spec)
#   * compute_on_startup=.false. on all FME AMs (warm-restart BFB)
#   * Non-FME AMs disabled to keep the tape minimal:
#       mpaso:  globalStats, timeSeriesStatsMonthly{,Min,Max}
#       mpassi: regionalStatistics, timeSeriesStatsDaily, timeSeriesStatsMonthly
#   * EAM fincl1 with FME-required fields and hist_file_storage_type='one_month'
#   * MPAS native streams gated to output_interval='none' (the .remapped.nc
#     files are the SamudrACE tape; native ne30 mesh files would be ~3-5x
#     more storage for the same data on a different grid)
#
# Leave user_nl() empty unless you need extra diagnostic output on top
# of the FME tape (e.g. a monthly fincl2 tape for sanity checks).

user_nl() {

# Example: add a monthly diagnostic tape alongside the FME tape.
#
# cat << 'EOF' >> user_nl_eam
# ! fincl2: monthly mean diagnostics (extra tape, supplements FME tape 1)
# fincl2 = 'PS','TS','PRECT','TMQ','FLUT','FSDS','FLDS'
# EOF

:
}

patch_mpas_streams() {

echo

}

######################################################
### Most users won't need to change anything below ###
######################################################

#-----------------------------------------------------
fetch_code() {

    if [ "${do_fetch_code,,}" != "true" ]; then
        echo $'\n----- Skipping fetch_code -----\n'
        return
    fi

    echo $'\n----- Starting fetch_code -----\n'
    local path=${CODE_ROOT}
    local repo=e3sm

    echo "Cloning $repo repository branch $BRANCH under $path"
    if [ -d "${path}" ]; then
        echo "ERROR: Directory already exists. Not overwriting"
        exit 20
    fi
    mkdir -p ${path}
    pushd ${path}

    git clone git@github.com:E3SM-Project/${repo}.git .

    rm -rf .git/hooks
    git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
    git config commit.template .git/hooks/commit.template

    git checkout ${BRANCH}

    if [ "${CHERRY}" != "" ]; then
        echo ----- WARNING: adding git cherry-pick -----
        for commit in "${CHERRY[@]}"
        do
            echo ${commit}
            git cherry-pick ${commit}
        done
        echo -------------------------------------------
    fi

    git submodule update --init --recursive

    popd
}

#-----------------------------------------------------
create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
        echo $'\n----- Skipping create_newcase -----\n'
        return
    fi

    echo $'\n----- Starting create_newcase -----\n'

    if [[ ${PELAYOUT} == custom-* ]]; then
        layout="M"
    else
        layout=${PELAYOUT}
    fi

    if [ ! -d "${FME_TESTMOD}" ]; then
        echo "ERROR: FME testmod not found at ${FME_TESTMOD}"
        echo "       Check CHECKOUT and CODE_ROOT settings."
        exit 30
    fi

    if [[ -z "$CASE_GROUP" ]]; then
        ${CODE_ROOT}/cime/scripts/create_newcase \
            --case ${CASE_NAME} \
            --output-root ${CASE_ROOT} \
            --script-root ${CASE_SCRIPTS_DIR} \
            --handle-preexisting-dirs u \
            --compset ${COMPSET} \
            --res ${RESOLUTION} \
            --machine ${MACHINE} \
            --project ${PROJECT} \
            --walltime ${WALLTIME} \
            --pecount ${layout} \
            --user-mods-dirs ${FME_TESTMOD}
    else
        ${CODE_ROOT}/cime/scripts/create_newcase \
            --case ${CASE_NAME} \
            --case-group ${CASE_GROUP} \
            --output-root ${CASE_ROOT} \
            --script-root ${CASE_SCRIPTS_DIR} \
            --handle-preexisting-dirs u \
            --compset ${COMPSET} \
            --res ${RESOLUTION} \
            --machine ${MACHINE} \
            --project ${PROJECT} \
            --walltime ${WALLTIME} \
            --pecount ${layout} \
            --user-mods-dirs ${FME_TESTMOD}
    fi

    if [ $? != 0 ]; then
      echo $'\nNote: if create_newcase failed because sub-directory already exists:'
      echo $'  * delete old case_script sub-directory'
      echo $'  * or set do_newcase=false\n'
      exit 35
    fi

}

#-----------------------------------------------------
case_setup() {

    if [ "${do_case_setup,,}" != "true" ]; then
        echo $'\n----- Skipping case_setup -----\n'
        return
    fi

    echo $'\n----- Starting case_setup -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    ./xmlchange EXEROOT=${CASE_BUILD_DIR}
    ./xmlchange RUNDIR=${CASE_RUN_DIR}

    ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING^^}
    ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

    if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then
      echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    user_nl

    ./case.setup --reset

    popd
}

#-----------------------------------------------------
custom_pelayout() {

if [[ ${PELAYOUT} == custom-* ]]; then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    if [ "${MACHINE}" == "pm-cpu" ]; then
        ncore=128
    elif [ "${MACHINE}" == "chrysalis" ]; then
        ncore=64
    elif [ "${MACHINE}" == "compy" ]; then
        ncore=40
    elif [ "${MACHINE}" == "anvil" ]; then
        ncore=36
    else
        echo 'ERROR: MACHINE = '${MACHINE}' is not supported for custom PE layout.'
        exit 400
    fi

    tmp=($(echo ${PELAYOUT} | tr "-" " "))
    nnodes=${tmp[1]}

    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange NTASKS=$(( $nnodes * $ncore ))
    ./xmlchange NTHRDS=1
    ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
    ./xmlchange MAX_TASKS_PER_NODE=$ncore
    popd

fi

}

#-----------------------------------------------------
case_build() {

    pushd ${CASE_SCRIPTS_DIR}

    if [ "${do_case_build,,}" != "true" ]; then

        echo $'\n----- case_build -----\n'

        if [ "${OLD_EXECUTABLE}" == "" ]; then
            if [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
                echo 'Skipping build because $do_case_build = '${do_case_build}
            else
                echo 'ERROR: $do_case_build = '${do_case_build}' but no executable exists for this case.'
                exit 297
            fi
        else
            if [ -x ${OLD_EXECUTABLE} ]; then
                echo 'Using $OLD_EXECUTABLE = '${OLD_EXECUTABLE}
                cp -fp ${OLD_EXECUTABLE} ${CASE_BUILD_DIR}/
            else
                echo 'ERROR: $OLD_EXECUTABLE = '$OLD_EXECUTABLE' does not exist or is not an executable file.'
                exit 297
            fi
        fi
        echo 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
        ./xmlchange BUILD_COMPLETE=TRUE

    else

        echo $'\n----- Starting case_build -----\n'

        if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
            ./xmlchange DEBUG=${DEBUG_COMPILE^^}
        fi

        ./case.build

        echo $'\n----- Preview namelists -----\n'
        ./preview_namelists

    fi

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    ./xmlchange RUN_STARTDATE=${START_DATE}
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}
    ./xmlchange BUDGETS=TRUE

    if (( RESUBMIT > 0 )); then
        ./xmlchange RESUBMIT=${RESUBMIT}
    fi

    if [ "${MODEL_START_TYPE,,}" == "initial" ]; then
        ./xmlchange RUN_TYPE="startup"
        ./xmlchange CONTINUE_RUN="FALSE"

    elif [ "${MODEL_START_TYPE,,}" == "continue" ]; then
        ./xmlchange CONTINUE_RUN="TRUE"

    elif [ "${MODEL_START_TYPE,,}" == "branch" ] || [ "${MODEL_START_TYPE,,}" == "hybrid" ]; then
        ./xmlchange RUN_TYPE=${MODEL_START_TYPE,,}
        ./xmlchange GET_REFCASE=${GET_REFCASE}
        ./xmlchange RUN_REFDIR=${RUN_REFDIR}
        ./xmlchange RUN_REFCASE=${RUN_REFCASE}
        ./xmlchange RUN_REFDATE=${RUN_REFDATE}
        echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE}
        echo '$RUN_REFDIR = '${RUN_REFDIR}
        echo '$RUN_REFCASE = '${RUN_REFCASE}
        echo '$RUN_REFDATE = '${START_DATE}

    else
        echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
        exit 380
    fi

    patch_mpas_streams

    popd
}

#-----------------------------------------------------
case_submit() {

    if [ "${do_case_submit,,}" != "true" ]; then
        echo $'\n----- Skipping case_submit -----\n'
        return
    fi

    echo $'\n----- Starting case_submit -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    ./case.submit

    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'

    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=$( basename -- "$0"; )
    local this_script_dir=$( dirname -- "$0"; )
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp "${this_script_dir}/${this_script_name}" ${script_provenance_dir}/${script_provenance_name}

}

#-----------------------------------------------------
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

#-----------------------------------------------------
main
