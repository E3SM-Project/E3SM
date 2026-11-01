#!/bin/bash -fe

# E3SM run script for FME (Full Model Emulation) production.
#
# Generates ACE/Samudra training data via the fme_output testmod, which
# configures online horizontal remapping, vertical coarsening, derived
# diagnostics, and column-integrated fields for EAM, MPAS-O, and MPAS-SI.
#
# Adapted from the standard E3SM run_e3sm.template.sh. Key differences:
#   * --user-mods-dirs added and points at cime_config/testmods_dirs/allactive/fme_output
#   * user_nl() is empty: the testmod is the single source of truth for
#     EAM/MPAS-O/MPAS-SI namelists -- including the FME AM enable/config and
#     the disable of non-FME AMs (globalStats, regionalStatistics,
#     timeSeriesStats*) that would otherwise add storage overhead. Add
#     fincl2/fincl3 here if you want extra diagnostic tapes alongside the
#     FME tape.
#   * Production restart cadence is 5 yr. With the 2026-05-01 restart fix
#     (append-mode reopen + accumulator sidecar + frame tracking +
#     compute_on_startup=.false.) the
#     model is BFB across restart, so the cadence choice is now about
#     wallclock budget and output-file size rather than restart-clobber
#     mitigation. 5-yr cadence with monthly file rotation keeps each
#     restart lookup small while still landing on a year-boundary, which
#     means leg N+1 always starts a fresh `*.YYYY-01.nc` (no append path
#     exercised in production).
#
#     See cime_config/testmods_dirs/allactive/fme_output/AGENTS.md for more information
#     and use it when modifying the FME code with an agent.
#   * The experiment is chosen by the FME_EXP preset (see the case statement
#     at the top of main). One script drives the whole family -- the F2010
#     PD/PI aerosol pair, F20TR, and the coupled WCYCL20TR / WCYCL1850 tapes --
#     so compset, start type, refcase, run length and PE layout live in one
#     table instead of being edited in place per run:
#         FME_EXP=F2010-PD ./run_e3sm.fme.sh
#         FME_EXP=WCYCL20TR ./run_e3sm.fme.sh

main() {

# --- Configuration flags ----
#
# Site-specific knobs (MACHINE, PROJECT, PELAYOUT, WALLTIME, RUN_REF*, and
# the scratch root via $PSCRATCH) are env-overridable. Defaults below are
# the pm-cpu/SamudrACE production combo. To run elsewhere, export e.g.
#   MACHINE=chrysalis PSCRATCH=/lcrc/group/e3sm/$USER/scratch ./run_e3sm.fme.sh
# COMPILER is intentionally omitted from create_newcase so each machine's
# default compiler is used; pin it explicitly here if you need a non-default.

# Default machine is ALCF crux (PBS, 128-core AMD Rome nodes, gnu). PROJECT
# is hardcoded rather than queried because crux is PBS, not Slurm -- the old
# `sacctmgr` lookup returns nothing here. crux's config_machines.xml already
# declares PROJECT=E3SMinput, so this just keeps the script self-consistent.
readonly MACHINE=${MACHINE:-crux}
readonly PROJECT="${PROJECT:-E3SMinput}"

# --- Experiment preset ------------------------------------------------------
# FME_EXP selects a named experiment. A preset sets only what actually differs
# between experiments -- compset, case name, start type, refcase, run length,
# PE layout, and an optional extra-namelist hook -- so adding a configuration
# means adding one branch to the case statement below, not editing the script.
# Everything a preset sets stays env-overridable, so the presets are defaults
# rather than a straitjacket.
#
#   F2010-PD    present-day aerosols, fixed-2010 forcing, data ocean (default)
#   F2010-PI    as F2010-PD, plus the pre-industrial aerosol overrides in
#               pi_aerosol_user_nl()
#   F20TR       transient atmosphere over prescribed SST / sea ice
#   WCYCL20TR   fully coupled historical, branched from a spun-up v3 LR case
#   WCYCL20TR-0301
#               coupled historical member 0301 (1940-2064, SSP2-4.5 after
#               2014), 0201's code/namelists/PE layout plus the 0301 output
#               spec; PILOT=TRUE runs only the 1940-1941 pilot leg
#   WCYCL1850   fully coupled piControl -- the original SamudrACE 100-yr tape
#
# Usage:
#   FME_EXP=F2010-PD ./run_e3sm.fme.sh
#   FME_EXP=WCYCL20TR ./run_e3sm.fme.sh
#   FME_EXP=WCYCL20TR-0301 PILOT=TRUE ./run_e3sm.fme.sh       # 2-yr pilot
#   FME_EXP=F2010-PI STOP_N=1 RESUBMIT=0 ./run_e3sm.fme.sh      # 1-yr shakeout
#   FME_EXP=F2010-PD STOP_OPTION=ndays STOP_N=10 RESUBMIT=0 ./run_e3sm.fme.sh
readonly FME_EXP="${FME_EXP:-F2010-PD}"
# PILOT=TRUE shortens a preset to its pilot leg (only WCYCL20TR-0301 has one).
readonly PILOT="${PILOT:-FALSE}"

# Preset outputs. Assigned plain here and frozen with readonly below, so an
# exported value always wins over the preset default.
exp_extra_user_nl=""   # shell function appended by user_nl(), or empty
exp_rest_n=""          # restart cadence in STOP_OPTION units
exp_refdir=""
exp_refcase=""
exp_refdate=""
exp_start_date=""
exp_walltime=""        # production walltime; empty = 24:00:00
exp_pe_tasks=""        # per-component "COMP:NTASKS:ROOTPE ..." on top of a custom-N layout
exp_xmlchange_run=""   # extra comma-separated env_run xmlchange settings, or empty

case "${FME_EXP^^}" in

  F2010-PD|F2010-PI)
    # F2010 = 2010_EAM%CMIP6_ELM%CNPRDCTCBCTOP_MPASSI%PRES_DOCN%DOM_MOSART_SGLC_SWAV:
    # prognostic EAM/ELM/MOSART over prescribed sea ice and a data ocean.
    # There is no MPAS-Ocean here, so the FME testmod's user_nl_mpaso blocks
    # land in a file CIME never reads and the MPAS-O FME streams are absent
    # from the tape; EAM, ELM, MOSART, MPAS-SI and coupler streams stay live.
    # OCN_NCPL defaults to $ATM_NCPL (48 on this grid) for a data ocean, so the
    # coupler-FME ocean-phase streams still flush on output boundaries and the
    # OCN_NCPL < ATM_NCPL regime of AGENTS.md gotcha #52 does not apply.
    #
    # Fixed-SST time slice: there is no ocean spinup to inherit and nothing to
    # branch from, so PD and PI both cold-start from identical F2010 initial
    # conditions and diverge only through the aerosol forcing.
    exp_scenario="${FME_EXP##*-}"
    exp_compset="F2010"
    exp_case_base="v3.LR.F2010.aigo.${exp_scenario^^}"
    exp_start_type="initial"
    exp_start_date="0001-01-01"
    # 10-yr legs x 3 submissions = 30 model years per scenario. Sized from
    # measured throughput rather than guessed: this configuration runs at
    # 17.8 s per model day (1.80 h per model year) on the layout below, so a
    # 10-yr leg is ~18.0 h of compute plus ~0.4 h of calendar-month file
    # rotation and ~0.1 h of annual restart writes -- ~18.6 h against crux's
    # 24 h workq-route cap, leaving ~5.4 h of margin. See AGENTS.md gotcha #65
    # for the cost model and for why a one-day timing run understates it by
    # 3.3x. Override with STOP_N / RESUBMIT.
    exp_stop_n="10"
    exp_resubmit="2"
    # Annual restarts: ~40 s and ~22.8 GB each, negligible against an 18.6 h
    # leg, and they give a recovery point every model year rather than only at
    # leg boundaries.
    exp_rest_n="1"
    # Reproduces the layout of the existing crux AMIP FME cases
    # (v3.LR.amip_*_aigoutput): 2880 ranks on all components, half-packed at
    # 64 ranks/node = 45 nodes. Half-packing is what those runs were tuned
    # with, so it is carried over rather than re-derived. 45 nodes also clears
    # crux's 9-node floor for the 24 h workq-route queue.
    exp_pelayout="custom-45-64"
    if [ "${exp_scenario^^}" == "PI" ]; then
        exp_extra_user_nl="pi_aerosol_user_nl"
    fi
    ;;

  F20TR)
    # Transient prescribed-SST atmosphere. Same component set as F2010, so the
    # same notes about MPAS-O absence and coupler streams apply.
    exp_compset="F20TR"
    exp_case_base="v3.LR.F20TR.aigo"
    exp_start_type="initial"
    exp_start_date="1850-01-01"
    exp_stop_n="5"
    exp_resubmit="1"
    exp_rest_n="5"                # one restart per leg (unvalidated preset)
    exp_pelayout="custom-45-64"   # as F2010, matching v3.LR.amip_*_aigoutput
    ;;

  WCYCL20TR)
    # Fully coupled historical, branched from a spun-up v3 LR historical run.
    # 'branch' = exact restart, bit-identical evolution with FME diagnostics
    # added on top; use MODEL_START_TYPE=hybrid if forcing/physics differ.
    exp_compset="WCYCL20TR"
    exp_case_base="v3.LR.historical_0201.aigo_cpl"
    exp_start_type="branch"
    exp_refdir="/pscratch/sd/m/mahf708/v3.LR.historical_0201/archive/rest/1940-01-01-00000"
    exp_refcase="v3.LR.historical_0201"
    exp_refdate="1940-01-01"
    exp_stop_n="5"
    exp_resubmit="4"
    exp_rest_n="5"                # one restart per leg (unvalidated preset)
    exp_pelayout="L"
    ;;

  WCYCL20TR-0301)
    # Coupled historical member 0301 for the ACE/Samudra training set: run
    # exactly like v3.LR.historical_0201.aigo_cpl (same code, namelists, PE
    # layout, coupling) except for the parent state and the 0301 output spec.
    # Any physics difference would make 0301 a different model.
    #
    # Parent: CMIP6 v3.LR.historical_0301 at 1940-01-01, staged by hand. Run
    #   0301_prep/scripts/prep_parent_restarts.sh <RUN_REFDIR>
    # on it first (total_water_ac fix applied to every earlier aigo parent,
    # rpointer and 5-day window alignment checks).
    #
    # Length: 1940-01-01 .. 2064-12-31 = 125 yr = 25 legs of 5 yr
    # (RESUBMIT=24). Legs end on Jan 1, which is always a 1850-anchored 5-day
    # window boundary, so no leg splits a window. eam.h0 for 2065-01 needs
    # nothing extra: the 2064-12-31 18Z..2065-01-01 00Z record is stamped at
    # the window end and lands in the last leg.
    exp_compset="WCYCL20TR"
    exp_case_base="v3.LR.historical_0301.aigo_cpl"
    exp_start_type="branch"
    exp_refdir="/eagle/E3SMinput/mahf708/scratch/crux/E3SMv3/v3.LR.historical_0301/archive/rest/1940-01-01-00000"
    exp_refcase="v3.LR.historical_0301"
    exp_refdate="1940-01-01"
    exp_stop_n="5"
    exp_resubmit="24"
    # 5-yearly Jan 1 sets as in 0201. The yearly May 1 / Nov 2 sets come from
    # restart_extra_mmdd in hist0301_user_nl(), not from REST_OPTION.
    exp_rest_n="5"
    exp_walltime="18:00:00"       # as 0201
    # 0201's layout, recovered from its env_mach_pes.xml: 30 full crux nodes,
    # pure MPI. ATM/CPL share nodes 0-21, LND/ROF 0-7 under them, ICE on 8-21,
    # OCN alone on 22-29. custom-30 gives 3840 tasks (also left on the stub
    # GLC/WAV/IAC/ESP components); exp_pe_tasks then places the active ones.
    # Equivalent to crux's tuned pesize="PM" entry in config_pesall.xml, but
    # spelled out so a PELAYOUT naming change can't silently move it.
    exp_pelayout="custom-30"
    exp_pe_tasks="CPL:2816:0 ATM:2816:0 LND:1024:0 ROF:1024:0 ICE:1792:1024 OCN:1024:2816"
    # 0201 coupling (30-min atm/lnd/ocn, 3-h river). These equal the defaults
    # for this compset/grid; pinned so a default change can't drift 0301 off
    # 0201 -- the coupler FME 6-h/1-d/5-d windows also assume them.
    exp_xmlchange_run="ATM_NCPL=48,LND_NCPL=48,OCN_NCPL=48,ROF_NCPL=8"
    exp_extra_user_nl="hist0301_user_nl"
    # Testmod knobs, read by fme_output/shell_commands during create_newcase,
    # so they must be in the environment before it runs. Still overridable.
    #   FME_EAM_PRES_LEVELS=FALSE  eam.h0 field list identical to 0101/0201
    #   FME_EAM_H1_ENABLE=TRUE     add the 5-day all-mean eam.h1 tape
    export FME_EAM_PRES_LEVELS="${FME_EAM_PRES_LEVELS:-FALSE}"
    export FME_EAM_H1_ENABLE="${FME_EAM_H1_ENABLE:-TRUE}"
    # PILOT=TRUE: only the 1940-1941 pilot (spec "Pilot validation"). It ends
    # 1942-01-01 (a window boundary) with an end-of-run restart set, written
    # because hist0301_user_nl() sets end_restart = .true. (the E3SM default is
    # .false.: REST_N=5 counts from the leg start, so a 2-yr leg otherwise ends
    # with no restart). To continue the same case afterwards, keep the
    # 5-yearly legs on 1945/1950/...: CONTINUE_RUN=TRUE, STOP_N=3 RESUBMIT=0
    # for 1942-1944, STOP_N=5 RESUBMIT=13 for 1945-2014, then the SSP2-4.5
    # switch (FME_0301_SSP245) and STOP_N=5 RESUBMIT=9 for 2015-2064.
    if [ "${PILOT^^}" == "TRUE" ]; then
        exp_stop_n="2"
        exp_resubmit="0"
    fi
    ;;

  WCYCL1850)
    # The original SamudrACE tape: 5-yr segments x 20 = 100 yr piControl,
    # branched from a spun-up piControl restart (a cold start would discard
    # ~50 yr of ocean spinup).
    exp_compset="WCYCL1850"
    exp_case_base="v3.LR.piControl.aigo"
    exp_start_type="branch"
    exp_refdir="/pscratch/sd/m/mahf708/v3.LR.piControl/archive/rest/0401-01-01-00000"
    exp_refcase="v3.LR.piControl"
    exp_refdate="0401-01-01"
    exp_stop_n="5"
    exp_resubmit="19"
    exp_rest_n="5"                # one restart per leg (unvalidated preset)
    exp_pelayout="L"
    ;;

  *)
    echo "ERROR: unknown FME_EXP '${FME_EXP}'."
    echo "       Known presets: F2010-PD, F2010-PI, F20TR, WCYCL20TR, WCYCL20TR-0301, WCYCL1850"
    echo "       Add a branch to the case statement in ${0##*/} for a new one."
    exit 5
    ;;
esac

readonly FME_EXTRA_USER_NL="${exp_extra_user_nl}"
# Per-component placement and extra env_run settings. '-' rather than ':-' so
# an explicit empty export (PE_TASKS= ...) switches the preset value off.
readonly PE_TASKS="${PE_TASKS-${exp_pe_tasks}}"
readonly EXTRA_XMLCHANGE_RUN="${EXTRA_XMLCHANGE_RUN-${exp_xmlchange_run}}"
readonly COMPSET="${COMPSET:-${exp_compset}}"
readonly RESOLUTION="${RESOLUTION:-ne30pg2_r05_IcoswISC30E3r5}"
# CASE_TAG is an optional dotted suffix (env-overridable) for launching
# several variants of one preset side-by-side without name collisions:
#   CASE_TAG=test1 FME_EXP=F2010-PD ./run_e3sm.fme.sh
#     -> v3.LR.F2010.aigo.PD.test1
readonly CASE_BASE="${CASE_BASE:-${exp_case_base}}"
readonly CASE_TAG="${CASE_TAG:-}"
readonly CASE_NAME="${CASE_BASE}${CASE_TAG:+.${CASE_TAG}}"
# readonly CASE_GROUP="samudrace_v3"

# Code and compilation
# BEFORE RUNNING: set CHECKOUT to a date string like 20260430
readonly CHECKOUT="latest"
readonly BRANCH="maint32/mahf708/fme/aigo"
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options, supplied by the preset above. 'initial' cold-starts; 'branch'
# and 'hybrid' additionally consume the RUN_REF* values, which the cold-start
# presets leave empty (they are inert on that path -- see runtime_options()).
readonly MODEL_START_TYPE="${MODEL_START_TYPE:-${exp_start_type}}"   # initial|continue|branch|hybrid
readonly RUN_REFDIR="${RUN_REFDIR:-${exp_refdir}}"
readonly RUN_REFCASE="${RUN_REFCASE:-${exp_refcase}}"
readonly RUN_REFDATE="${RUN_REFDATE:-${exp_refdate}}"
# A branch run starts where its refcase stopped; a cold start uses the preset's
# own calendar origin.
readonly START_DATE="${START_DATE:-${exp_start_date:-${RUN_REFDATE}}}"

# GET_REFCASE=TRUE asks CIME to copy refcase restart files from RUN_REFDIR.
# Set FALSE if you pre-staged them into CASE_RUN_DIR yourself. Ignored unless
# MODEL_START_TYPE is branch or hybrid.
readonly GET_REFCASE="${GET_REFCASE:-TRUE}"

# Set paths.
# PSCRATCH is a NERSC-ism; on crux the equivalent is the machine's
# CIME_OUTPUT_ROOT, /eagle/$PROJECT/$USER/scratch/crux. Export PSCRATCH to
# override (e.g. when running this script at NERSC).
readonly SCRATCH_ROOT="${PSCRATCH:-/eagle/${PROJECT}/${USER}/scratch/crux}"
readonly CASE_ROOT="${SCRATCH_ROOT}/E3SMv3/${CASE_NAME}"
# CODE_ROOT points at the shared crux checkout already on
# maint32/mahf708/fme/aigo, so every experiment compiles from one tree. Cases
# that share a compset can also share an executable via OLD_EXECUTABLE.
# do_fetch_code is false to match; set it true and point CODE_ROOT at a fresh
# path to get the old self-contained per-case clone instead.
readonly CODE_ROOT="${CODE_ROOT:-/lus/eagle/projects/E3SMinput/mahf708/scratch/crux/mahf708-e3sm}"

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

  # Production: length and PE layout come from the FME_EXP preset, so each
  # experiment carries its own tape size (10 yr for the F2010 aerosol pair,
  # 25 yr for the historical branch, 100 yr for the piControl tape). Override
  # either in the environment, e.g.
  #   FME_EXP=F2010-PI STOP_N=1 RESUBMIT=0 ./run_e3sm.fme.sh
  # Legs land on year boundaries, so each leg opens a fresh monthly file
  # (`*.YYYY-01.remapped.nc`) and the mid-window append path stays unexercised.
  # The append-mode + sidecar machinery (AGENTS.md #29) keeps mid-window
  # restarts safe regardless; the cadence just avoids exercising it.
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  # crux queue policy (config_batch.xml) is strict: `debug` is 1-8 nodes capped
  # at 2 h, `workq-route` is 9-184 nodes capped at 24 h, so a production leg
  # must ask for at least 9 nodes.
  #
  # The F presets pin an explicit custom layout because crux's tuned entry for
  # this grid in config_pesall.xml matches only MPASO-bearing (WCYCL) compsets
  # -- an F compset would otherwise fall back to a generic default. The WCYCL
  # presets keep the named 'L' layout, which does hit crux's tuned entry.
  readonly PELAYOUT=${PELAYOUT:-${exp_pelayout}}
  readonly WALLTIME=${WALLTIME:-${exp_walltime:-24:00:00}}
  # Cadence knobs are overridable like everything else the preset supplies.
  # STOP_OPTION was previously pinned to 'nyears' while STOP_N was
  # overridable, which made short shakeout legs (STOP_OPTION=ndays STOP_N=10)
  # impossible through this script -- you had to drive CIME directly.
  # REST_OPTION follows STOP_OPTION so a restart cadence is expressed in the
  # same units as the leg unless explicitly overridden.
  readonly STOP_OPTION="${STOP_OPTION:-nyears}"
  readonly STOP_N="${STOP_N:-${exp_stop_n}}"
  readonly REST_OPTION="${REST_OPTION:-${STOP_OPTION}}"
  readonly REST_N="${REST_N:-${exp_rest_n}}"
  readonly RESUBMIT="${RESUBMIT:-${exp_resubmit}}"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history
readonly HIST_OPTION="nyears"
readonly HIST_N="5"

# Reuse an existing executable only when the component set matches. The old
# WCYCL exe has MPAS-Ocean compiled in and F2010 swaps that for a data ocean,
# so the F presets must build from scratch. Within a preset family it is safe
# and worthwhile: once F2010-PD has built, set
#   OLD_EXECUTABLE=<PD case>/build/e3sm.exe
# for F2010-PI -- same source, same compset, only namelists differ -- and skip
# the second compile.
readonly OLD_EXECUTABLE="${OLD_EXECUTABLE:-}"

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=${DO_CASE_SUBMIT:-false}

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
    # A missing parent set (the 0301 parent is staged by hand) is fatal only
    # when this invocation submits; create + build may legitimately run before
    # staging, and the WCYCL20TR/WCYCL1850 presets keep their NERSC refdirs.
    if [ "${GET_REFCASE^^}" == "TRUE" ] && [ ! -d "${RUN_REFDIR}" ]; then
        if [ "${do_case_submit,,}" == "true" ]; then
            echo "ERROR: RUN_REFDIR=${RUN_REFDIR} does not exist."
            echo "       Stage the ${RUN_REFCASE} ${RUN_REFDATE} restart set there first."
            exit 11
        fi
        echo "WARNING: RUN_REFDIR=${RUN_REFDIR} does not exist yet."
        echo "         Stage the ${RUN_REFCASE} ${RUN_REFDATE} restart set before case.submit."
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
#
# Per-experiment namelist differences do NOT belong here -- they go in a
# dedicated hook function named by the preset's exp_extra_user_nl, which
# user_nl() invokes below. That keeps "what this experiment is" next to the
# preset instead of scattered behind if-statements.

user_nl() {

# Example: add a monthly diagnostic tape alongside the FME tape.
#
# cat << 'EOF' >> user_nl_eam
# ! fincl2: monthly mean diagnostics (extra tape, supplements FME tape 1)
# fincl2 = 'PS','TS','PRECT','TMQ','FLUT','FSDS','FLDS'
# EOF

# Preset-supplied namelist hook (empty for most presets).
if [ -n "${FME_EXTRA_USER_NL}" ]; then
    echo "Applying experiment namelist hook: ${FME_EXTRA_USER_NL}"
    ${FME_EXTRA_USER_NL}
fi

:
}

# ---------------------------------------------------------------------------
# Pre-industrial aerosol overrides for the F2010-PI half of the aerosol pair.
# ---------------------------------------------------------------------------
#
# F2010-PD is the stock F2010 compset with no namelist overrides at all, so
# everything that makes PI a pre-industrial-aerosol run lives in this one
# function. It runs from inside CASE_SCRIPTS_DIR, before `case.setup --reset`,
# so it can append to user_nl_eam (and any other user_nl_* file) directly.
#
# NOT YET FILLED IN. The guard below makes F2010-PI fail loudly at setup rather
# than quietly producing a second copy of PD -- two identical runs labelled PD
# and PI is the one failure mode that would silently poison the forcing
# difference. Replace the guard with the aerosol settings, e.g. the 1850
# emission streams and prescribed-oxidant files:
#
#   cat << 'EOF' >> user_nl_eam
#   ! Pre-industrial (1850) aerosol and precursor emissions
#   ! ext_frc_specifier  = ...
#   ! srf_emis_specifier = ...
#   ! tracer_cnst_file   = ...
#   EOF
#
pi_aerosol_user_nl() {

    # Ask CIME for the machine's inputdata root instead of hardcoding one.
    # The upstream spec these lists came from was written against Chrysalis
    # (/lcrc/group/e3sm/data/inputdata) and every path had to be retranslated
    # by hand for crux (/grand/E3SMinput/data) -- doing it this way means the
    # next machine costs nothing.
    local din emis
    din=$(./xmlquery --value DIN_LOC_ROOT)
    emis="${din}/atm/cam/chem/trop_mozart_aero/emis"

    # Aerosol-ONLY perturbation: GHG concentrations, prescribed SST/sea ice,
    # ozone, solar and land use all stay at their F2010 values. Only the
    # aerosol and aerosol-precursor emissions move to 1850, so the PD-PI
    # difference is the aerosol effective radiative forcing and nothing else.
    #
    # Two file families are in play, which is expected rather than an
    # inconsistency:
    #   * chem_gases + DMS + E90 + aircraft NO2 use purpose-built
    #     "2010-as-1850" climatologies (2010 seasonal cycle, 1850 magnitudes).
    #   * DECK_ne30 aerosol species use the transient 1850-2014 files, from
    #     which *_cycle_yr = 1850 selects the 1850 slice.
    # Both are driven CYCLICAL at cycle_yr 1850.
    #
    # The species lists below are identical to what F2010-PD generates
    # (22 surface, 10 elevated) -- the pair must differ in emission VALUES
    # only. Verify with:
    #   diff <(sed -n "/srf_emis_specifier/,/srf_emis_type/p" PD/CaseDocs/atm_in) \
    #        <(sed -n "/srf_emis_specifier/,/srf_emis_type/p" PI/CaseDocs/atm_in)

    cat << EOF >> user_nl_eam
!--- Pre-industrial (1850) aerosol and precursor emissions (F2010-PI)

ext_frc_cycle_yr  = 1850
ext_frc_specifier = 'NO2    -> ${emis}/chem_gases/2degrees/emissions-cmip6_NO2_aircraft_vertical_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'SO2         -> ${emis}/DECK_ne30/cmip6_mam4_so2_elev_1850-2014_c180205.nc',
        'SOAG0       -> ${emis}/DECK_ne30/emissions-cmip6_e3sm_SOAG0_elev_1850-2014_1.9x2.5_c20230201.nc',
        'bc_a4       -> ${emis}/DECK_ne30/cmip6_mam4_bc_a4_elev_1850-2014_c180205.nc',
        'num_a1      -> ${emis}/DECK_ne30/cmip6_mam4_num_a1_elev_1850-2014_c180205.nc',
        'num_a2      -> ${emis}/DECK_ne30/cmip6_mam4_num_a2_elev_1850-2014_c180205.nc',
        'num_a4      -> ${emis}/DECK_ne30/cmip6_mam4_num_a4_elev_1850-2014_c180205.nc',
        'pom_a4      -> ${emis}/DECK_ne30/cmip6_mam4_pom_a4_elev_1850-2014_c180205.nc',
        'so4_a1      -> ${emis}/DECK_ne30/cmip6_mam4_so4_a1_elev_1850-2014_c180205.nc',
        'so4_a2      -> ${emis}/DECK_ne30/cmip6_mam4_so4_a2_elev_1850-2014_c180205.nc'
ext_frc_type      = 'CYCLICAL'

srf_emis_cycle_yr  = 1850
srf_emis_specifier = 'C2H4     -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_C2H4_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'C2H6     -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_C2H6_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'C3H8     -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_C3H8_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'CH2O     -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_CH2O_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'CH3CHO   -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_CH3CHO_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'CH3COCH3 -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_CH3COCH3_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'CO       -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_CO_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'ISOP     -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_ISOP_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'ISOP_VBS -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_ISOP_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'C10H16   -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_MTERP_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'SOAG0    -> ${emis}/DECK_ne30/emissions-cmip6_e3sm_SOAG0_surf_1850-2014_1.9x2.5_c20230201.nc',
        'NO       -> ${emis}/chem_gases/2degrees/emissions-cmip6_e3sm_NO_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc',
        'DMS      -> ${emis}/DMSflux.2010-as-1850.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20190220.nc',
        'SO2      -> ${emis}/DECK_ne30/cmip6_mam4_so2_surf_1850-2014_c180205.nc',
        'bc_a4    -> ${emis}/DECK_ne30/cmip6_mam4_bc_a4_surf_1850-2014_c180205.nc',
        'num_a1   -> ${emis}/DECK_ne30/cmip6_mam4_num_a1_surf_1850-2014_c180205.nc',
        'num_a2   -> ${emis}/DECK_ne30/cmip6_mam4_num_a2_surf_1850-2014_c180205.nc',
        'num_a4   -> ${emis}/DECK_ne30/cmip6_mam4_num_a4_surf_1850-2014_c180205.nc',
        'pom_a4   -> ${emis}/DECK_ne30/cmip6_mam4_pom_a4_surf_1850-2014_c180205.nc',
        'so4_a1   -> ${emis}/DECK_ne30/cmip6_mam4_so4_a1_surf_1850-2014_c180205.nc',
        'so4_a2   -> ${emis}/DECK_ne30/cmip6_mam4_so4_a2_surf_1850-2014_c180205.nc',
        'E90      -> ${emis}/chem_gases/2degrees/emissions_E90_surface_2010-as-1850_clim_1.9x2.5_c20230213.nc'
srf_emis_type      = 'CYCLICAL'
EOF

}

# ---------------------------------------------------------------------------
# Namelist additions for WCYCL20TR-0301 (coupled historical member 0301).
# ---------------------------------------------------------------------------
#
# Runs from CASE_SCRIPTS_DIR after create_newcase, i.e. AFTER the fme_output
# testmod has written its user_nl_* blocks, so anything here is appended
# below them. CIME's user_nl parser lets a later scalar setting win over an
# earlier one, so a repeat of a testmod variable here overrides it, not errors.
#
# 1. SSP2-4.5 forcing from 2015 (EAM and ELM), ONLY when FME_0301_SSP245=TRUE.
#    Copied verbatim from the hand edits made to 0201's user_nl_eam /
#    user_nl_elm on Aug 8 (recovered in 0301_prep/recovered_0201/), which
#    0201 applied when its run reached 2015: 32 of the 38 files start in 2015
#    (GHG 2015-2500, emissions 2015-2100, land use 2015-2100), so EAM aborts
#    at a 1940 start with them ("chem_surfvals_set_all: data time index is
#    out of bounds"). 1940-2014 therefore runs on the compset's historical
#    defaults, and at the 2015-01-01 leg boundary the block is appended to
#    the case's user_nl_eam / user_nl_elm. Either re-run this function with
#    FME_0301_SSP245=TRUE from the case dir, or append the saved copies
#    case_scripts/ssp245_from2015.user_nl_{eam,elm}. The $DIN_LOC_ROOT tokens
#    are left for build-namelist to expand (quoted heredocs), as 0201 had them.
#
# 2. Coupler / driver output switches for the 0301 spec (user_nl_cpl feeds
#    both drv_in and the cpl_fme_inparm group). Kept here rather than in the
#    testmod so every other experiment keeps its current coupler behaviour.
#      restart_extra_mmdd       extra full restart sets at 00Z on May 1 and
#                               Nov 2 every year, on top of REST_OPTION's
#                               5-yearly Jan 1 (both dates are 1850-anchored
#                               5-day window boundaries; Nov 1 is not)
#      cpl_fme_state_mean       1-day/5-day coupler state fields written as
#                               window means instead of end-of-window instants
#      cpl_fme_inst_enable      add the *5Dinst (and *1Dinst) end-of-window
#                               state snapshot streams
#      cpl_fme_window_anchor_ymd  count the coupler FME windows from 1850-01-01
#                               explicitly instead of from the branch date
hist0301_user_nl() {

    if [ "${FME_0301_SSP245^^}" == "TRUE" ]; then
    cat << 'EOF' >> user_nl_eam

!--- WCYCL20TR-0301: CMIP6 SSP2-4.5 forcing after 2014 (verbatim from 0201.aigo_cpl)
  chlorine_loading_file		= '$DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/Linoz_Chlorine_Loading_CMIP6_Hist_SSP245_0003-2503_c20200808.nc'
  ext_frc_specifier		= 'NO2         -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_NO2_aircraft_vertical_2015-2100_1.9x2.5_c20240219.nc',
          'SO2         -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_so2_volc_elev_2015-2100_c240331.nc',
          'SOAG0       -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_SOAG0_elev_2015-2100_1.9x2.5_c20240219.nc',
          'bc_a4       -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_bc_a4_elev_2015-2100_c200716.nc',
          'num_a1      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_num_a1_elev_2015-2100_c200716.nc',
          'num_a2      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_num_a2_elev_2015-2100_c200716.nc',
          'num_a4      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_num_a4_elev_2015-2100_c200716.nc',
          'pom_a4      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_pom_a4_elev_2015-2100_c200716.nc',
          'so4_a1      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_so4_a1_elev_2015-2100_c200716.nc',
          'so4_a2      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_so4_a2_elev_2015-2100_c200716.nc'
  linoz_data_file		= 'linv3_1849-2101_CMIP6_Hist_SSP245_10deg_58km_c20230705.nc'
  srf_emis_specifier		= 'C10H16 -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_MTERP_surface_2015-2100_1.9x2.5_c20240219.nc',
          'C2H4      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_C2H4_surface_2015-2100_1.9x2.5_c20240219.nc',
          'C2H6      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_C2H6_surface_2015-2100_1.9x2.5_c20240219.nc',
          'C3H8      -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_C3H8_surface_2015-2100_1.9x2.5_c20240219.nc',
          'CH2O   -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_CH2O_surface_2015-2100_1.9x2.5_c20240219.nc',
          'CH3CHO    -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_CH3CHO_surface_2015-2100_1.9x2.5_c20240219.nc',
          'CH3COCH3  -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_CH3COCH3_surface_2015-2100_1.9x2.5_c20240219.nc',
          'CO     -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_CO_surface_2015-2100_1.9x2.5_c20240219.nc',
          'DMS    -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc',
          'E90       -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/emissions_E90_surface_1750-2101_1.9x2.5_c20231222.nc',
          'ISOP   -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_ISOP_surface_2015-2100_1.9x2.5_c20240219.nc',
          'ISOP_VBS -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_ISOP_surface_2015-2100_1.9x2.5_c20240219.nc',
          'NO     -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_NO_surface_2015-2100_1.9x2.5_c20240219.nc',
          'SO2    -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_so2_surf_2015-2100_c200716.nc',
          'SOAG0  -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/emissions-cmip6_ssp245_e3sm_SOAG0_surf_2015-2100_1.9x2.5_c20240219.nc',
          'bc_a4  -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_bc_a4_surf_2015-2100_c200716.nc',
          'num_a1 -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_num_a1_surf_2015-2100_c200716.nc',
          'num_a2 -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_num_a2_surf_2015-2100_c200716.nc',
          'num_a4 -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_num_a4_surf_2015-2100_c200716.nc',
          'pom_a4 -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_pom_a4_surf_2015-2100_c200716.nc',
          'so4_a1 -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_so4_a1_surf_2015-2100_c200716.nc',
          'so4_a2 -> $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP245_ne30/cmip6_ssp245_mam4_so4_a2_surf_2015-2100_c200716.nc'
  tracer_cnst_file		= 'oxid_SSP245_1.9x2.5_L70_1849-2101_c20240228.nc'
  bndtvghg		= '$DIN_LOC_ROOT/atm/cam/ggas/GHG_CMIP_SSP245-1-2-1_Annual_Global_2015-2500_c20200807.nc'
EOF

    cat << 'EOF' >> user_nl_elm

!--- WCYCL20TR-0301: CMIP6 SSP2-4.5 land use, N deposition and population density (verbatim from 0201.aigo_cpl)
 flanduse_timeseries = '${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_ssp2_rcp45_simyr2015-2100_c240408.nc'
 stream_fldfilename_ndep = '${DIN_LOC_ROOT}/lnd/clm2/ndepdata/fndep_elm_cbgc_exp_simyr1849-2101_1.9x2.5_ssp245_c240903.nc'
 stream_year_last_ndep = 2100
 stream_fldfilename_popdens = '${DIN_LOC_ROOT}/lnd/clm2/firedata/elmforc.Li_20181205_mod_hist_SSP2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c240906.nc'
 stream_year_last_popdens = 2100
EOF

    fi
    if [ "${FME_0301_SSP245_ONLY^^}" == "TRUE" ]; then
        return   # 2015 switch on an existing case: SSP block only
    fi

    cat << 'EOF' >> user_nl_cpl

!--- WCYCL20TR-0301: yearly May 1 / Nov 2 restarts and window-mean coupler streams
restart_extra_mmdd        = 501, 1102
! E3SM's end_restart default is .false.; without it a leg whose length is not
! a multiple of REST_N (the 2-yr pilot) ends with no restart set.
end_restart               = .true.
cpl_fme_state_mean        = .true.
cpl_fme_inst_enable       = .true.
cpl_fme_window_anchor_ymd = 18500101
EOF

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

    # Preset-pinned env_run settings (e.g. the 0301 coupling frequencies).
    # Set before user_nl/case.setup so preview_namelists already reflects them.
    if [ -n "${EXTRA_XMLCHANGE_RUN}" ]; then
        echo "Applying preset env_run settings: ${EXTRA_XMLCHANGE_RUN}"
        ./xmlchange "${EXTRA_XMLCHANGE_RUN}"
    fi

    # COSP simulator intentionally NOT enabled. The FME/SamudrACE tape carries
    # no COSP fields, so '-cosp' (docosp=.true.) would run the simulator every
    # few radiation steps for 100 years and produce nothing on tape. Re-add
    # "--id CAM_CONFIG_OPTS --append --val='-cosp'" here if COSP output is ever
    # wanted.

    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    user_nl

    ./case.setup --reset

    popd
}

#-----------------------------------------------------
custom_pelayout() {

if [[ ${PELAYOUT} == custom-* ]]; then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    if [ "${MACHINE}" == "crux" ]; then
        ncore=128
    elif [ "${MACHINE}" == "pm-cpu" ]; then
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
    # Optional third field sets ranks-per-node, e.g. custom-45-64 = 45 nodes
    # at 64 ranks/node (half-packed). Omit it to pack the machine's full core
    # count, i.e. custom-12 on crux == custom-12-128.
    nranks=${tmp[2]:-$ncore}

    echo "  ${nnodes} nodes x ${nranks} ranks/node = $(( $nnodes * $nranks )) tasks"

    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange NTASKS=$(( $nnodes * $nranks ))
    ./xmlchange NTHRDS=1
    ./xmlchange MAX_MPITASKS_PER_NODE=$nranks
    ./xmlchange MAX_TASKS_PER_NODE=$nranks

    # Optional per-component placement on top of the uniform layout above,
    # "COMP:NTASKS:ROOTPE" entries (the WCYCL20TR-0301 preset uses it to
    # reproduce 0201). Applied only with the preset's own PELAYOUT: a list
    # sized for 30 nodes would be wrong under any other node count.
    if [ -n "${PE_TASKS}" ] && [ "${PELAYOUT}" == "${exp_pelayout}" ]; then
        local entry comp ntasks rootpe
        # The uniform block above leaves ROOTPE as --pecount M's fallback entry
        # set it; start every component (stubs included) from 0.
        ./xmlchange ROOTPE=0
        for entry in ${PE_TASKS}; do
            IFS=: read -r comp ntasks rootpe <<< "${entry}"
            if (( rootpe + ntasks > nnodes * nranks )); then
                echo "ERROR: PE_TASKS entry ${entry} exceeds ${nnodes}x${nranks} tasks"
                exit 401
            fi
            echo "  ${comp}: ${ntasks} tasks at rootpe ${rootpe}"
            ./xmlchange NTASKS_${comp}=${ntasks},ROOTPE_${comp}=${rootpe},NTHRDS_${comp}=1
        done
    fi
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
