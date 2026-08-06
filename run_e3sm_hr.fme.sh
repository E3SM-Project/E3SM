#!/bin/bash -fe

# E3SM run script for FME (Full Model Emulation) production, layered on top of
# the v3.HR 1950 historical respinup.
#
# This script BUILDS ON TOP OF ac.jwolfe's
#   /lcrc/group/e3sm2/ac.jwolfe/run.v3.HR.1950_historical_respinup_0051.sh
# It hybrid-restarts from the same spun-up ocean/ice/land state (refcase
# 20260702.v3.HR.1950-CMIP7_respinpup_0041_elm0001 @ 0051-01-01), runs the same
# WCYCL20TR-CMIP7 transient configuration with the same physics and forcing, and
# produces the SAME diagnostic output as the respinup -- then ADDS the FME
# (ACE/Samudra) training tapes on top, non-destructively.
#
# How the two output sets are kept from colliding:
#   * EAM: the respinup owns history tapes 1-6 (cosp/aerosol/chem). The FME
#     field list is relocated to tape 7 via `export FME_EAM_TAPE=7`, so the
#     fme_output testmod writes fincl7 / horiz_remap_file(7) while this script's
#     user_nl() owns empty_htapes and the 7-element per-tape arrays. Disjoint
#     tape indices => no clobber (ptapes=17 allows it).
#   * MPAS-O / MPAS-SI: the FME analysis members write their own remapped
#     streams; this script RE-ENABLES the non-FME AMs (globalStats,
#     timeSeriesStats*, highFrequencyOutput, ...) that the testmod disables for
#     storage, restoring the respinup's ocean/ice output alongside the FME one.
#   * ELM / MOSART: the testmod sets hist_empty_htapes=.true.; this script
#     overrides back to .false. so the respinup's land/river tapes survive.
#   * cpl: the FME coupler-native streams (cpl_fme_*) don't touch any respinup
#     coupler output.
#
# Operational knobs (machine, PE layout, restart cadence, resubmit, submit
# toggle, case naming) are LEFT as in the standard FME run script. The hybrid
# refcase + compset + full science namelists are brought over from the respinup.

main() {

# --- Configuration flags ----
#
# Site-specific knobs (MACHINE, PROJECT, PELAYOUT, WALLTIME, RUN_REF*, and
# the scratch root via $PSCRATCH) are env-overridable.
# NOTE: the hybrid refcase below lives on chrysalis (/lcrc/group/e3sm2/...),
# so run this on chrysalis, e.g.
#   MACHINE=chrysalis PELAYOUT=custom-225 ./run_e3sm_hr.fme.sh
# (or export RUN_REFDIR/RUN_REFCASE to a copy staged on your machine).

readonly MACHINE=${MACHINE:-pm-cpu}
readonly PROJECT="${PROJECT:-$(sacctmgr show user $USER format=DefaultAccount | tail -n1 | tr -d ' ')}"

readonly COMPSET="WCYCL20TR-CMIP7"
readonly RESOLUTION="ne120pg2_r025_RRSwISC6to18E3r5"
# BEFORE RUNNING: set CASE_BASE to a descriptive identifier. CASE_TAG is
# an optional dotted suffix (env-overridable) for launching several
# variants side-by-side without name collisions.
readonly CASE_BASE="${CASE_BASE:-v3.HR.testprod.aigo}"
readonly CASE_TAG="${CASE_TAG:-}"
readonly CASE_NAME="${CASE_BASE}${CASE_TAG:+.${CASE_TAG}}"
# readonly CASE_GROUP="samudrace_v3"

# Code and compilation
# BEFORE RUNNING: set CHECKOUT to a date string like 20260430
readonly CHECKOUT="latest"
readonly BRANCH="maint32/mahf708/fme/aigo"
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options.
#
# Hybrid-restart from the v3.HR 1950 historical respinup (brought over from
# ac.jwolfe's run.v3.HR.1950_historical_respinup_0051.sh). This reuses ~50 yr
# of spun-up ocean/ice/land state; the atmosphere reinitializes with the
# WCYCL20TR-CMIP7 transient forcing. START_DATE maps refcase model-year 0051
# to calendar 1950 and the run marches forward historically from there.
readonly MODEL_START_TYPE="hybrid"   # 'initial', 'continue', 'branch', 'hybrid'

# Hybrid refcase (from the respinup). All are env-overridable.
readonly RUN_REFDIR="${RUN_REFDIR:-/lcrc/group/e3sm2/ac.jwolfe/E3SMv3_dev/20260702.v3.HR.1950-CMIP7_respinpup_0041_elm0001/archive/rest/0051-01-01-00000}"
readonly RUN_REFCASE="${RUN_REFCASE:-20260702.v3.HR.1950-CMIP7_respinpup_0041_elm0001}"
readonly RUN_REFDATE="${RUN_REFDATE:-0051-01-01}"
readonly START_DATE="${START_DATE:-1950-01-01}"

# GET_REFCASE=TRUE asks CIME to copy refcase restart files from RUN_REFDIR.
# Set to FALSE if you've pre-staged the files into CASE_RUN_DIR yourself.
readonly GET_REFCASE=TRUE

# Set paths
readonly CASE_ROOT="${PSCRATCH}/E3SMv3/${CASE_NAME}"
# CODE_ROOT defaults to a self-contained checkout under CASE_ROOT so each
# case is self-describing. Override by exporting CODE_ROOT=/path/to/checkout.
readonly CODE_ROOT="${CODE_ROOT:-${CASE_ROOT}/code/${CHECKOUT}}"

readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# FME testmod directory (applied via --user-mods-dirs at create_newcase time)
readonly FME_TESTMOD="${CODE_ROOT}/cime_config/testmods_dirs/allactive/fme_output"

# Put the FME EAM field list on history tape 7 so it does not collide with the
# respinup's tapes 1-6. The fme_output shell_commands reads FME_EAM_TAPE and
# emits fincl7 + horiz_remap_file(7); this script's user_nl() owns empty_htapes
# and the 7-element per-tape arrays (see user_nl() below). Keep the tape-7 slot
# in those arrays in sync with the FME_EAM_* defaults (6-hourly, one_month,
# mfilt 1500, avgflag 'I').
export FME_EAM_TAPE=7

# Define type of run
#  short tests: 'XS_2x5_ndays', 'S_1x10_ndays', etc. (same scheme as upstream)
#  or 'production' for the full FME + respinup production tape
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

  # Production run. Cadence/PE/walltime are the FME script's own operational
  # settings (left unchanged); only the science config is inherited from the
  # respinup.
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT=${PELAYOUT:-L}
  readonly WALLTIME=${WALLTIME:-24:00:00}
  readonly STOP_OPTION="nyears"
  readonly STOP_N="20"
  readonly REST_OPTION="nyears"
  readonly REST_N="5"
  readonly RESUBMIT="4"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history
readonly HIST_OPTION="nyears"
readonly HIST_N="5"

readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=false
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
# The fme_output testmod (applied via --user-mods-dirs) runs FIRST at
# create_newcase time and appends the FME namelist to user_nl_*. This user_nl()
# runs AFTER (in case_setup) and appends the respinup's science namelists on
# top. Because namelist is last-wins for scalars and element-wise for arrays,
# the two coexist cleanly:
#   * EAM: FME writes fincl7 / horiz_remap_file(7) / vcoarsen+derived+aerocom
#     globals (empty_htapes and the per-tape arrays are delegated to us because
#     FME_EAM_TAPE=7 != 1). We provide empty_htapes and the 7-element arrays
#     (respinup tapes 1-6 + FME tape 7).
#   * MPAS-O/SI: FME disables the non-FME AMs; we RE-ENABLE them here to restore
#     the respinup's ocean/ice output.
#   * ELM/MOSART: FME sets *_empty_htapes=.true.; we override back to .false.

user_nl() {

# ---------------------------------------------------------------------------
# EAM: respinup tapes 1-6 + physics/forcing. FME field list is tape 7 (testmod).
# ---------------------------------------------------------------------------
cat << EOF >> user_nl_eam

 cosp_lite = .true.

 empty_htapes = .true.

 ! 7 history tapes: respinup tapes 1-6 (below) + FME field tape 7 (written by
 ! the fme_output testmod via FME_EAM_TAPE=7). The tape-7 slot in each array
 ! MUST match the FME_EAM_* env defaults (FME_EAM_OUTPUT_HOURS=6 -> nhtfrq -6;
 ! FME_EAM_MFILT=1500; FME_EAM_STORAGE=one_month; FME_EAM_AVGFLAG='I'). The
 ! testmod owns fincl7, horiz_remap_file(7), and the vcoarsen/derived/aerocom
 ! globals -- do NOT set those here.
 avgflag_pertape = 'A','A','I','A','I','I','I'
 nhtfrq = 0,-24,-6,-3,-1,0,-6
 mfilt  = 1,30,120,240,720,1,1500
 hist_file_storage_type = 'num_snapshots','num_snapshots','num_snapshots','num_snapshots','num_snapshots','num_snapshots','one_month'

 fincl1 = 'AODALL','AODBC','AODDUST','AODPOM','AODSO4','AODSOA','AODSS','AODVIS',
          'CLDLOW','CLDMED','CLDHGH','CLDTOT',
          'CLDHGH_CAL','CLDLOW_CAL','CLDMED_CAL','CLD_MISR','CLDTOT_CAL',
          'CLMODIS','FISCCP1_COSP','FLDS','FLDSC','FLNS','FLNSC','FLNT','FLUT',
          'FLUTC','FSDS','FSDSC','FSNS','FSNSC','FSNT','FSNTOA','FSNTOAC',
          'ICEFRAC','LANDFRAC','LWCF','OCNFRAC','OMEGA','PRECC','PRECL','PRECSC','PRECSL','PS','PSL','Q',
          'QFLX','QREFHT','RELHUM','SCO','SHFLX','SOLIN','SWCF','T','TAUX','TAUY','TCO',
          'TGCLDLWP','TMQ','TREFHT','TREFMNAV','TREFMXAV','TS','U','U10','V','Z3',
          'O3','LHFLX',
          'ABURDENSO4_STR','ABURDENSO4_TRO',
          'ABURDENSO4','ABURDENBC','ABURDENDUST','ABURDENMOM','ABURDENPOM','ABURDENSEASALT',
          'ABURDENSOA','AODSO4_STR',
          'H2OLNZ',
          'PHIS','CLOUD','TGCLDIWP','TGCLDCWP','AREL',
 fincl2 = 'PS', 'FLUT','PRECT','PRECC','U200','V200','U850','V850',
 fincl4 = 'PRECT'
 fincl6 = 'CO_2DMSD','NO2_2DMSD','NO_2DMSD','O3_2DMSD','O3_2DMSD_trop'

 ! -- chemUCI settings ------------------
 history_chemdyg_summary = .true.
 history_gaschmbudget_2D = .false.
 history_gaschmbudget_2D_levels = .false.
 history_gaschmbudget_num = 6 !! no impact if  history_gaschmbudget_2D = .false.

 ! -- MAM5 settings ---------------------
 is_output_interactive_volc = .true.

 clubb_c8       =  5.2
 nucleate_ice_subgrid   = 1.40

EOF

# ---------------------------------------------------------------------------
# ELM: restore land output (override the testmod's hist_empty_htapes=.true.)
# plus respinup physics/forcing.
# ---------------------------------------------------------------------------
cat << EOF >> user_nl_elm

 ! Override the fme_output testmod's hist_empty_htapes=.true. so the respinup
 ! land history tapes below are produced.
 hist_empty_htapes = .false.

 hist_dov2xy = .true.,.true.
 hist_fexcl1 = 'AGWDNPP','ALTMAX_LASTYEAR','AVAIL_RETRANSP','AVAILC','BAF_CROP',
               'BAF_PEATF','BIOCHEM_PMIN_TO_PLANT','CH4_SURF_AERE_SAT','CH4_SURF_AERE_UNSAT','CH4_SURF_DIFF_SAT',
               'CH4_SURF_DIFF_UNSAT','CH4_SURF_EBUL_SAT','CH4_SURF_EBUL_UNSAT','CMASS_BALANCE_ERROR','cn_scalar',
               'COL_PTRUNC','CONC_CH4_SAT','CONC_CH4_UNSAT','CONC_O2_SAT','CONC_O2_UNSAT',
               'cp_scalar','CWDC_HR','CWDC_LOSS','CWDC_TO_LITR2C','CWDC_TO_LITR3C',
               'CWDC_vr','CWDN_TO_LITR2N','CWDN_TO_LITR3N','CWDN_vr','CWDP_TO_LITR2P',
               'CWDP_TO_LITR3P','CWDP_vr','DWT_CONV_CFLUX_DRIBBLED','F_CO2_SOIL','F_CO2_SOIL_vr',
               'F_DENIT_vr','F_N2O_DENIT','F_N2O_NIT','F_NIT_vr','FCH4_DFSAT',
               'FINUNDATED_LAG','FPI_P_vr','FPI_vr','FROOTC_LOSS','HR_vr',
               'LABILEP_TO_SECONDP','LABILEP_vr','LAND_UPTAKE','LEAF_MR','leaf_npimbalance',
               'LEAFC_LOSS','LEAFC_TO_LITTER','LFC2','LITR1_HR','LITR1C_TO_SOIL1C',
               'LITR1C_vr','LITR1N_TNDNCY_VERT_TRANS','LITR1N_TO_SOIL1N','LITR1N_vr','LITR1P_TNDNCY_VERT_TRANS',
               'LITR1P_TO_SOIL1P','LITR1P_vr','LITR2_HR','LITR2C_TO_SOIL2C','LITR2C_vr',
               'LITR2N_TNDNCY_VERT_TRANS','LITR2N_TO_SOIL2N','LITR2N_vr','LITR2P_TNDNCY_VERT_TRANS','LITR2P_TO_SOIL2P',
               'LITR2P_vr','LITR3_HR','LITR3C_TO_SOIL3C','LITR3C_vr','LITR3N_TNDNCY_VERT_TRANS',
               'LITR3N_TO_SOIL3N','LITR3N_vr','LITR3P_TNDNCY_VERT_TRANS','LITR3P_TO_SOIL3P','LITR3P_vr',
               'M_LITR1C_TO_LEACHING','M_LITR2C_TO_LEACHING','M_LITR3C_TO_LEACHING','M_SOIL1C_TO_LEACHING','M_SOIL2C_TO_LEACHING',
               'M_SOIL3C_TO_LEACHING','M_SOIL4C_TO_LEACHING','NDEPLOY','NEM','nlim_m',
               'o2_decomp_depth_unsat','OCCLP_vr','PDEPLOY','PLANT_CALLOC','PLANT_NDEMAND',
               'PLANT_NDEMAND_COL','PLANT_PALLOC','PLANT_PDEMAND','PLANT_PDEMAND_COL','plim_m',
               'POT_F_DENIT','POT_F_NIT','POTENTIAL_IMMOB','POTENTIAL_IMMOB_P','PRIMP_TO_LABILEP',
               'PRIMP_vr','PROD1P_LOSS','QOVER_LAG','RETRANSN_TO_NPOOL','RETRANSP_TO_PPOOL',
               'SCALARAVG_vr','SECONDP_TO_LABILEP','SECONDP_TO_OCCLP','SECONDP_vr','SMIN_NH4_vr',
               'SMIN_NO3_vr','SMINN_TO_SOIL1N_L1','SMINN_TO_SOIL2N_L2','SMINN_TO_SOIL2N_S1','SMINN_TO_SOIL3N_L3',
               'SMINN_TO_SOIL3N_S2','SMINN_TO_SOIL4N_S3','SMINP_TO_SOIL1P_L1','SMINP_TO_SOIL2P_L2','SMINP_TO_SOIL2P_S1',
               'SMINP_TO_SOIL3P_L3','SMINP_TO_SOIL3P_S2','SMINP_TO_SOIL4P_S3','SMINP_vr','SOIL1_HR',
               'SOIL1C_TO_SOIL2C','SOIL1C_vr','SOIL1N_TNDNCY_VERT_TRANS','SOIL1N_TO_SOIL2N','SOIL1N_vr',
               'SOIL1P_TNDNCY_VERT_TRANS','SOIL1P_TO_SOIL2P','SOIL1P_vr','SOIL2_HR','SOIL2C_TO_SOIL3C',
               'SOIL2C_vr','SOIL2N_TNDNCY_VERT_TRANS','SOIL2N_TO_SOIL3N','SOIL2N_vr','SOIL2P_TNDNCY_VERT_TRANS',
               'SOIL2P_TO_SOIL3P','SOIL2P_vr','SOIL3_HR','SOIL3C_TO_SOIL4C','SOIL3C_vr',
               'SOIL3N_TNDNCY_VERT_TRANS','SOIL3N_TO_SOIL4N','SOIL3N_vr','SOIL3P_TNDNCY_VERT_TRANS','SOIL3P_TO_SOIL4P',
               'SOIL3P_vr','SOIL4_HR','SOIL4C_vr','SOIL4N_TNDNCY_VERT_TRANS','SOIL4N_TO_SMINN',
               'SOIL4N_vr','SOIL4P_TNDNCY_VERT_TRANS','SOIL4P_TO_SMINP','SOIL4P_vr','SOLUTIONP_vr',
               'TCS_MONTH_BEGIN','TCS_MONTH_END','TOTCOLCH4','water_scalar','WF',
               'wlim_m','WOODC_LOSS','WTGQ'
 hist_fincl1 = 'SNOWDP','COL_FIRE_CLOSS','NPOOL','PPOOL','TOTPRODC'
 hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA'
 hist_mfilt = 1,30
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'
 check_finidat_year_consistency = .false.
 check_finidat_fsurdat_consistency = .false.

 !finidat = '/lcrc/group/e3sm/data/inputdata/lnd/clm2/initdata_map/elmi.CNPRDCTCBCTOP.r025_RRSwISC6to18E3r5.1950-01-01-00000.c20260412.nc'
 check_finidat_pct_consistency = .false.

 model_year_align_ndep = 1849
 stream_fldfilename_ndep = '/lcrc/group/e3sm/data/inputdata/lnd/clm2/ndepdata/ndep_from_CMIP7_184901-202312_v2.0_1.9x2.5_aave_c260515.nc'
 stream_year_first_ndep = 1849
 stream_year_last_ndep = 2023

 fsurdat = '/lcrc/group/e3sm/data/inputdata/lnd/clm2/surfdata_map/surfdata_r025_simyr1950_cmip7_c260404.nc'
 flanduse_timeseries = '/lcrc/group/e3sm/data/inputdata/lnd/clm2/surfdata_map/landuse.timeseries_r025_hist_simyr1950-2024_cmip7_c260404.nc'
EOF

# ---------------------------------------------------------------------------
# MOSART: restore river output (override the testmod's rtmhist_empty_htapes).
# ---------------------------------------------------------------------------
cat << EOF >> user_nl_mosart

 rtmhist_empty_htapes = .false.
EOF

# ---------------------------------------------------------------------------
# MPAS-Ocean: respinup physics + RE-ENABLE the non-FME analysis members that
# the fme_output testmod disables, so the respinup's ocean output is restored
# alongside the FME remapped streams.
# ---------------------------------------------------------------------------
cat << EOF >> user_nl_mpaso

 ! Re-enable non-FME analysis members disabled by the fme_output testmod so the
 ! respinup's ocean diagnostic output is produced (FME streams are separate).
 config_AM_globalStats_enable               = .true.
 config_AM_timeSeriesStatsMonthly_enable    = .true.
 config_AM_timeSeriesStatsMonthlyMin_enable = .true.
 config_AM_timeSeriesStatsMonthlyMax_enable = .true.
 config_AM_highFrequencyOutput_enable       = .true.
 config_AM_meridionalHeatTransport_enable   = .true.
 config_AM_oceanHeatContent_enable          = .true.

 ! for efficiency
 config_compute_active_tracer_budgets = .false.

 config_cvmix_kpp_langmuir_mixing_opt = 'NONE'
 config_cvmix_kpp_langmuir_entrainment_opt = 'NONE'
 config_cvmix_background_scheme = 'latitude-dependent'
 config_cvmix_background_diffusion = 3.0e-5
 config_cvmix_background_diffusion_max_latitude = -50

EOF

# ---------------------------------------------------------------------------
# MPAS-SeaIce: RE-ENABLE regionalStatistics + timeSeriesStatsMonthly (disabled
# by the testmod), keep daily off for the spinup (respinup setting).
# ---------------------------------------------------------------------------
cat << EOF >> user_nl_mpassi

 ! Re-enable non-FME sea-ice AMs disabled by the fme_output testmod.
 config_AM_regionalStatistics_enable     = .true.
 config_AM_timeSeriesStatsMonthly_enable = .true.

 ! turn off daily output for the spinup
 config_am_timeseriesstatsdaily_enable = .false.

EOF

}

patch_mpas_streams() {

echo
echo 'Modifying MPAS streams files'
pushd ../run

# change streams.ocean file
patch streams.ocean << EOF
--- streams.ocean
+++ streams.ocean
@@ -331,7 +331,7 @@
         filename_interval="00-01-00_00:00:00"
         reference_time="01-01-01_00:00:00"
         output_interval="00-00-05_00:00:00"
-        clobber_mode="append"
+        clobber_mode="truncate"
         packages="highFrequencyOutputAMPKG">

     <var name="xtime"/>
@@ -526,1 +526,1 @@
-        output_interval="00-00-01_00:00:00"
+        output_interval="01-00-00_00:00:00"
@@ -661,7 +661,6 @@
     <var name="windStressZonal"/>
     <var name="windStressMeridional"/>
     <var name="ekeCorrection"/>
-    <var name="frazilLayerThicknessTendency"/>
     <var_array name="avgValueWithinOceanRegion"/>
     <var_array name="avgValueWithinOceanLayerRegion"/>
     <var_array name="avgValueWithinOceanVolumeRegion"/>
@@ -691,7 +690,6 @@
     <var name="longWaveHeatFluxDown"/>
     <var name="seaIceHeatFlux"/>
     <var name="shortWaveHeatFlux"/>
-    <var name="frazilTemperatureTendency"/>
     <var name="evaporationFlux"/>
     <var name="seaIceSalinityFlux"/>
     <var name="seaIceFreshWaterFlux"/>
@@ -701,23 +699,13 @@
     <var name="rainFlux"/>
     <var name="snowFlux"/>
     <var name="bottomLayerShortwaveTemperatureFlux"/>
-    <var_array name="activeTracersTend"/>
     <var name="salinitySurfaceRestoringTendency"/>
-    <var_array name="activeTracerHorizontalAdvectionTendency"/>
-    <var_array name="activeTracerVerticalAdvectionTendency"/>
-    <var_array name="activeTracerVertMixTendency"/>
-    <var_array name="activeTracerHorMixTendency"/>
-    <var_array name="activeTracerSurfaceFluxTendency"/>
-    <var_array name="temperatureShortWaveTendency"/>
-    <var_array name="activeTracerNonLocalTendency"/>
     <var name="areaCellGlobal"/>
     <var name="areaEdgeGlobal"/>
     <var name="areaTriangleGlobal"/>
     <var name="volumeCellGlobal"/>
     <var name="volumeEdgeGlobal"/>
     <var name="CFLNumberGlobal"/>
-    <var name="vertDiffTopOfCell"/>
-    <var name="vertViscTopOfCell"/>
     <var name="boundaryLayerDepth"/>
     <var name="columnIntegratedSpeed"/>
     <var name="landIceFreshwaterFlux"/>
@@ -726,25 +714,16 @@
     <var name="landIceHeatFlux"/>
     <var name="heatFluxToLandIce"/>
     <var name="landIceFrictionVelocity"/>
-    <var name="velocityTidalRMS"/>
     <var_array name="landIceBoundaryLayerTracers"/>
     <var_array name="landIceInterfaceTracers"/>
     <var name="mocStreamvalLatAndDepth"/>
     <var name="mocStreamvalLatAndDepthRegion"/>
     <var name="binBoundaryMocStreamfunction"/>
     <var name="surfaceBuoyancyForcing"/>
-    <var name="tendLayerThickness"/>
     <var name="pressureAdjustedSSH"/>
     <var name="SSHSquared"/>
     <var name="velocityZonalSquared"/>
     <var name="velocityMeridionalSquared"/>
-    <var name="normalVelocitySquared"/>
-    <var name="velocityZonalTimesTemperature"/>
-    <var name="velocityMeridionalTimesTemperature"/>
-    <var name="normalVelocityTimesTemperature"/>
-    <var name="velocityZonalTimesSalinity"/>
-    <var name="velocityMeridionalTimesSalinity"/>
-    <var name="normalVelocityTimesSalinity"/>
     <var name="oceanHeatContentSfcToBot"/>
     <var name="oceanHeatContentSfcTo700m"/>
     <var name="oceanHeatContent700mTo2000m"/>
@@ -756,13 +735,9 @@
     <var name="activeTracerHorAdvectionMLTend"/>
     <var name="activeTracerVertMixMLTend"/>
     <var name="activeTracersML"/>
-    <var name="BruntVaisalaFreqTop"/>
     <var name="bruntVaisalaFreqML"/>
     <var name="activeTracersTendML"/>
-    <var_array name="activeTracerVerticalAdvectionTopFlux"/>
-    <var_array name="activeTracerHorizontalAdvectionEdgeFlux"/>
     <var_array name="totalFreshWaterTemperatureFlux"/>
-    <var name="velocityTidalRMS"/>
 </stream>

 <stream name="timeSeriesStatsMonthlyMaxOutput"
EOF

# copy to SourceMods
cp streams.ocean ../case_scripts/SourceMods/src.mpaso/

popd

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

    # Enable the COSP simulator (the respinup's fincl1 includes COSP cloud
    # diagnostics: FISCCP1_COSP, CLD*_CAL, CLD_MISR, CLMODIS -- and cosp_lite is
    # set in user_nl_eam). Skipped for a data atmosphere.
    if [ `./xmlquery --value COMP_ATM` == "datm" ]; then
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
        ./xmlchange CONTINUE_RUN="FALSE"
        ./xmlchange GET_REFCASE=${GET_REFCASE}
        ./xmlchange RUN_REFDIR=${RUN_REFDIR}
        ./xmlchange RUN_REFCASE=${RUN_REFCASE}
        ./xmlchange RUN_REFDATE=${RUN_REFDATE}
        echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE}
        echo '$RUN_REFDIR = '${RUN_REFDIR}
        echo '$RUN_REFCASE = '${RUN_REFCASE}
        echo '$RUN_REFDATE = '${RUN_REFDATE}
        # mpas streams files aren't updated with branch/hybrid info until
        # preview_namelists runs; regenerate them before patch_mpas_streams so
        # the patch applies to the final streams.ocean.
        echo $'\n----- Preview namelists -----\n'
        ./preview_namelists

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
