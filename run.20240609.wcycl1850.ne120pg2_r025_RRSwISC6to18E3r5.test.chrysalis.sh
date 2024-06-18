#!/bin/bash -fe

# E3SM Coupled Model Group run_e3sm script template.
#
# Bash coding style inspired by:
# http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

main() {

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=chrysalis
readonly PROJECT="e3sm"

# Simulation
readonly COMPSET="WCYCL1850"
readonly RESOLUTION="ne120pg2_r025_RRSwISC6to18E3r5"
readonly NL_MAPS=true   ### nonlinear maps for tri-grid, waiting for trfvnp2 maps
readonly CASE_NAME="20240609.piCtl.ne120pg2_r025_RRSwISC6to18E3r5.chrysalis.test1" 

# Code and compilation
readonly CHECKOUT="20240604-rrswisc6to18e3r5-elmpatch"
readonly BRANCH="master" #based on master 93e511d57c5619d0888d63be4a1e99bda64cf8c9 
readonly CHERRY=()
readonly DEBUG_COMPILE=false

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="0001-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=FALSE
readonly RUN_REFDIR="/lcrc/group/e3sm2/ac.xzheng/E3SMv3_dev/20231117.v3b02.piControl.chrysalis/archive/rest/0101-01-01-00000"
readonly RUN_REFCASE="20231117.v3b02.piControl.chrysalis"  
readonly RUN_REFDATE="0101-01-01"

# Set paths
readonly CODE_ROOT="/lcrc/group/e3sm2/ac.xzheng/E3SMv3_dev/code/${CHECKOUT}"
readonly CASE_ROOT="/lcrc/group/e3sm2/ac.xzheng/E3SMv3_dev/${CASE_NAME}"

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_1x10_ndays', 'XS_2x5_ndays', 'S_1x10_ndays', 'M_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation


#readonly run='custom-96-production'
#readonly run='custom-56_1x10_ndays'
#readonly run='custom-22_1x31_ndays'
#readonly run='custom-220_1x1_nmonths'
#readonly run='custom-88_1x1_nmonths'
readonly run='production'
#readonly run='custom-44_1x1_nmonths'
#readonly run='custom-110_1x1_nmonths'
#readonly run='custom-220_1x1_nmonths'

if [[ "${run}" != *"production"* ]]; then
  echo "setting up Short test simulations: ${run}"
  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}

  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/tests/${run}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/tests/${run}/run
  readonly PELAYOUT=${layout}
  readonly WALLTIME="03:00:00"
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false

else
  echo "setting up ${run}"
  # Production simulation
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="custom-80"
  readonly WALLTIME="48:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="2"
  readonly REST_OPTION="nmonths"
  readonly REST_N="3"
  readonly RESUBMIT="0"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history 
readonly HIST_OPTION="nyears"
readonly HIST_N="1"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=false

# --- Now, do the work ---

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

# Custom PE layout
custom_pelayout

# Setup
case_setup

# Build
case_build

# Configure runtime options
runtime_options

# Copy script into case_script directory for provenance
copy_script

# Submit
case_submit

# All done
echo $'\n----- All done -----\n'

}

# =======================
# Custom user_nl settings
# =======================

user_nl() {

cat << EOF >> user_nl_eam
 cosp_lite = .true.

 empty_htapes = .true.

 avgflag_pertape = 'A','A','A','A','I','I'
 nhtfrq = 0,-24,-6,-3,-1,0
 mfilt  = 1,30,120,240,720,1

 fincl1 = 'AODALL','AODBC','AODDUST','AODPOM','AODSO4','AODSOA','AODSS','AODVIS',
          'CLDLOW','CLDMED','CLDHGH','CLDTOT',
          'CLDHGH_CAL','CLDLOW_CAL','CLDMED_CAL','CLD_MISR','CLDTOT_CAL',
          'CLMODIS','FISCCP1_COSP','FLDS','FLNS','FLNSC','FLNT','FLUT',
          'FLUTC','FSDS','FSDSC','FSNS','FSNSC','FSNT','FSNTOA','FSNTOAC','FSNTC',
          'ICEFRAC','LANDFRAC','LWCF','OCNFRAC','OMEGA','PRECC','PRECL','PRECSC','PRECSL','PS','PSL','Q',
          'QFLX','QREFHT','RELHUM','SCO','SHFLX','SOLIN','SWCF','T','TAUX','TAUY','TCO',
          'TGCLDLWP','TMQ','TREFHT','TREFMNAV','TREFMXAV','TS','U','U10','V','Z3',
          'dst_a1DDF','dst_a3DDF','dst_c1DDF','dst_c3DDF','dst_a1SFWET','dst_a3SFWET','dst_c1SFWET','dst_c3SFWET',
          'O3','LHFLX',
          'O3_2DTDA_trop','O3_2DTDB_trop','O3_2DTDD_trop','O3_2DTDE_trop','O3_2DTDI_trop','O3_2DTDL_trop',
          'O3_2DTDN_trop','O3_2DTDO_trop','O3_2DTDS_trop','O3_2DTDU_trop','O3_2DTRE_trop','O3_2DTRI_trop',
          'O3_SRF','NO_2DTDS','NO_TDLgt','NO2_2DTDD','NO2_2DTDS','NO2_TDAcf','CO_SRF','TROPE3D_P','TROP_P',
          'CDNUMC','SFDMS','so4_a1_sfgaex1','so4_a2_sfgaex1','so4_a3_sfgaex1','so4_a5_sfgaex1','soa_a1_sfgaex1',
          'soa_a2_sfgaex1','soa_a3_sfgaex1','GS_soa_a1','GS_soa_a2','GS_soa_a3','AQSO4_H2O2','AQSO4_O3',
          'SFSO2','SO2_CLXF','SO2','DF_SO2','AQ_SO2','GS_SO2','WD_SO2','ABURDENSO4_STR','ABURDENSO4_TRO',
          'ABURDENSO4','ABURDENBC','ABURDENDUST','ABURDENMOM','ABURDENPOM','ABURDENSEASALT',
          'ABURDENSOA','AODSO4_STR','AODSO4_TRO',
          'EXTINCT','AODABS','AODABSBC','CLDICE','CLDLIQ','CLD_CAL_TMPLIQ','CLD_CAL_TMPICE','Mass_bc_srf',
          'Mass_dst_srf','Mass_mom_srf','Mass_ncl_srf','Mass_pom_srf','Mass_so4_srf','Mass_soa_srf','Mass_bc_850',
          'Mass_dst_850','Mass_mom_850','Mass_ncl_850','Mass_pom_850','Mass_so4_850','Mass_soa_850','Mass_bc_500',
          'Mass_dst_500','Mass_mom_500','Mass_ncl_500','Mass_pom_500','Mass_so4_500','Mass_soa_500','Mass_bc_330',
          'Mass_dst_330','Mass_mom_330','Mass_ncl_330','Mass_pom_330','Mass_so4_330','Mass_soa_330','Mass_bc_200',
          'Mass_dst_200','Mass_mom_200','Mass_ncl_200','Mass_pom_200','Mass_so4_200','Mass_soa_200',
          'O3_2DTDD','O3_2DCIP','O3_2DCIL','CO_2DTDS','CO_2DTDD','CO_2DCEP','CO_2DCEL','NO_2DTDD',
          'FLNTC','SAODVIS',
          'H2OLNZ',
          'dst_a1SF','dst_a3SF',
          'PHIS','CLOUD','TGCLDIWP','TGCLDCWP','AREL',
          'CLDTOT_ISCCP','MEANCLDALB_ISCCP','MEANPTOP_ISCCP','CLD_CAL',
          'CLDTOT_CAL_LIQ','CLDTOT_CAL_ICE','CLDTOT_CAL_UN',
          'CLDHGH_CAL_LIQ','CLDHGH_CAL_ICE','CLDHGH_CAL_UN',
          'CLDMED_CAL_LIQ','CLDMED_CAL_ICE','CLDMED_CAL_UN',
          'CLDLOW_CAL_LIQ','CLDLOW_CAL_ICE','CLDLOW_CAL_UN',
          'CLWMODIS','CLIMODIS'

 fincl2 = 'PS', 'FLUT','PRECT','U200','V200','U850','V850',
          'TCO','SCO','TREFHTMN','TREFHTMX','TREFHT','QREFHT'
 fincl3 = 'PS', 'PSL','PRECT','TUQ','TVQ','UBOT','VBOT','TREFHT','FLUT','OMEGA500','TBOT','U850','V850','U200','V200','T200','T500','Z700'
 fincl4 = 'PRECT'
 fincl5 = 'O3_SRF'
 fincl6 = 'CO_2DMSD','NO2_2DMSD','NO_2DMSD','O3_2DMSD','O3_2DMSD_trop'

 ! -- chemUCI settings ------------------
 history_chemdyg_summary = .true.
 history_gaschmbudget_2D = .false.
 history_gaschmbudget_2D_levels = .false.
 history_gaschmbudget_num = 6 !! no impact if  history_gaschmbudget_2D = .false.

 ! -- MAM5 settings ------------------    
 is_output_interactive_volc = .true.                                                                     

 ! Parameter changes for HR

 cld_macmic_num_steps           =  3

 ! Turn mountain stress on

  do_tms = .true.
                                                                                                               
EOF

cat << EOF >> user_nl_elm
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
              'SMINP_TO_SOIL3P_L3','SMINP_TO_SOIL3P_S2','SMINP_TO_SOIL4P_S3','SMINP_vr','SOIL1_HR','SOIL1C_TO_SOIL2C','SOIL1C_vr','SOIL1N_TNDNCY_VERT_TRANS','SOIL1N_TO_SOIL2N','SOIL1N_vr',
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
 hist_mfilt = 1,365
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'

 check_finidat_year_consistency = .false.
 check_finidat_fsurdat_consistency = .false.

 ! if finidat from a different period is specified
 ! check_finidat_pct_consistency   = .false.

 !--- land BGC spin-up initial conditions ---, pending
 finidat='/lcrc/group/e3sm/data/inputdata/lnd/clm2/initdata_map/elmi.CNPRDCTCBCTOP.r025_RRSwISC6to18E3r5.1850.nc'
EOF

cat << EOF >> user_nl_mpaso
 config_btr_dt = '0000_00:00:05'
 config_mom_del2 = 100.0
 config_use_mom_del2 = .true.
EOF
}

# =====================================
# Customize MPAS stream files if needed
# =====================================
patch_mpas_streams() {

echo
echo 'Modifying MPAS streams files' 
pushd ../run
# change streams.ocean file
patch streams.ocean << EOF
--- streams.ocean       2024-06-06 13:40:29.317030813 -0500
+++ streams.ocean       2024-06-06 13:41:25.340746391 -0500
@@ -9,7 +9,7 @@
                   type="input"
                   io_type="pnetcdf,cdf5"
                   input_interval="initial_only"
-                  filename_template="/lcrc/group/e3sm/data/inputdata/ocn/mpas-o/RRSwISC6to18E3r5/mpaso.RRSwISC6to18E3r5.20240327.nc"
+                  filename_template="/lcrc/group/e3sm/data/inputdata/ocn/mpas-o/RRSwISC6to18E3r5/mpaso.RRSwISC6to18E3r5.rstFromG-chrysalis.20240603.nc"
 />

 <!--
EOF
# change streams.seaice file
patch streams.seaice << EOF
--- streams.seaice    # Original file content
+++ streams.seaice    # Modified file content
@@ -8,7 +8,7 @@
                   type="input"
                   io_type="pnetcdf,cdf5"
                   filename_interval="none"
-                  filename_template="/lcrc/group/e3sm/data/inputdata/ice/mpas-seaice/RRSwISC6to18E3r5/mpassi.RRSwISC6to18E3r5.20240327.nc"
+                  filename_template="/lcrc/group/e3sm/data/inputdata/ice/mpas-seaice/RRSwISC6to18E3r5/mpassi.RRSwISC6to18E3r5.rstFromG-chrysalis.20240603.nc"
                   input_interval="initial_only" />

 <!--
@@ -35,7 +35,7 @@
 <immutable_stream name="restart_ic"
                   type="input"
                   io_type="pnetcdf,cdf5"
-                  filename_template="/lcrc/group/e3sm/data/inputdata/ice/mpas-seaice/RRSwISC6to18E3r5/mpassi.RRSwISC6to18E3r5.20240327.nc"
+                  filename_template="/lcrc/group/e3sm/data/inputdata/ice/mpas-seaice/RRSwISC6to18E3r5/mpassi.RRSwISC6to18E3r5.rstFromG-chrysalis.20240603.nc"
                   filename_interval="none"
                   input_interval="initial_only" />

EOF

# copy to SourceMods
cp streams.ocean ../case_scripts/SourceMods/src.mpaso/
cp streams.seaice ../case_scripts/SourceMods/src.mpassi/

popd

}
# =====================================================
# Custom PE layout: custom-N where N is number of nodes
# =====================================================

custom_pelayout(){

if [[ ${PELAYOUT} == custom-* ]];
then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    # Number of cores per node (machine specific)
    if [ "${MACHINE}" == "chrysalis" ]; then
        ncore=64
        hthrd=2  # hyper-threading
    elif [ "${MACHINE}" == "pm-cpu" ]; then
        ncore=128
        hthrd=1  # hyper-threading
    else
        echo 'ERROR: MACHINE = '${MACHINE}' is not supported for current custom PE layout setting.'
        exit 400
    fi

    # Extract number of nodes
    tmp=($(echo ${PELAYOUT} | tr "-" " "))
    nnodes=${tmp[1]}

    # Applicable to all custom layouts
    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange NTASKS=1
    ./xmlchange NTHRDS=1
    ./xmlchange ROOTPE=0
    ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
    ./xmlchange MAX_TASKS_PER_NODE=$(( $ncore * $hthrd))

    # Layout-specific customization
    if [ "${nnodes}" == "100" ]; then

       echo Using custom 100 nodes layout

       ### Current defaults for L
      ./xmlchange CPL_NTASKS=5440
      ./xmlchange ATM_NTASKS=5440
      ./xmlchange OCN_NTASKS=960
      ./xmlchange OCN_ROOTPE=5440

      ### Added by Xue for tri-grid
        ./xmlchange LND_NTASKS=1600
        ./xmlchange ROF_NTASKS=1600
        ./xmlchange ICE_NTASKS=3840
        ./xmlchange LND_ROOTPE=3840
        ./xmlchange ROF_ROOTPE=3840
    elif [ "${nnodes}" == "56" ]; then
       echo Using 56 nodes stacked layout on Chrysalis
      ./xmlchange CPL_NTASKS=3584
      ./xmlchange ATM_NTASKS=3584
      ./xmlchange OCN_NTASKS=3584
      ./xmlchange LND_NTASKS=3584
      ./xmlchange ROF_NTASKS=3584
      ./xmlchange ICE_NTASKS=3584

    elif [ "${nnodes}" == "56F" ]; then
       echo Using 56 nodes layout on Chrysalis

      ./xmlchange CPL_NTASKS=3168
      ./xmlchange ATM_NTASKS=3168
      ./xmlchange OCN_NTASKS=416
      ./xmlchange OCN_ROOTPE=3168 

      ### Added by Xue for tri-grid
        ./xmlchange LND_NTASKS=1120  
        ./xmlchange ROF_NTASKS=1120 
        ./xmlchange ICE_NTASKS=2048 
        ./xmlchange LND_ROOTPE=2048 
        ./xmlchange ROF_ROOTPE=2048 
    elif [ "${nnodes}" == "20" ]; then
       echo Using  20 nodes stocked layout on Chrysalis
      ./xmlchange CPL_NTASKS=1280
      ./xmlchange ATM_NTASKS=1280
      ./xmlchange OCN_NTASKS=1280
      ./xmlchange LND_NTASKS=1280
      ./xmlchange ROF_NTASKS=1280
      ./xmlchange ICE_NTASKS=1280

    elif [ "${nnodes}" == "96" ]; then

       echo Using custom 96 nodes layout

      ./xmlchange CPL_NTASKS=9216
      ./xmlchange ATM_NTASKS=9216
      ./xmlchange OCN_NTASKS=3072
      ./xmlchange OCN_ROOTPE=9216

      ### Added by Xue for tri-grid
        ./xmlchange LND_NTASKS=3072
        ./xmlchange ROF_NTASKS=3072
        ./xmlchange ICE_NTASKS=6144
        ./xmlchange LND_ROOTPE=6144
        ./xmlchange ROF_ROOTPE=6144

     elif [ "${nnodes}" == "11" ]; then

      ### Current defaults for M layout on chrysalis for F-Cases
      # The custom-22 was previously tested for bi-grid with L80-QBO

      ./xmlchange CPL_NTASKS=1408
      ./xmlchange ATM_NTASKS=1408   #let atm run on all CPU
      ./xmlchange OCN_NTASKS=192
      ./xmlchange OCN_ROOTPE=1216

      ./xmlchange LND_NTASKS=1152
      ./xmlchange ROF_NTASKS=1152
      ./xmlchange ICE_NTASKS=256
      ./xmlchange LND_ROOTPE=256
      ./xmlchange ROF_ROOTPE=256

     elif [ "${nnodes}" == "22old" ]; then

      ### Current defaults for M layout on chrysalis for F-Cases
      # The custom-22 was previously tested for bi-grid with L80-QBO

      ./xmlchange CPL_NTASKS=1408
      ./xmlchange ATM_NTASKS=1408   #let atm run on all CPU
      ./xmlchange OCN_NTASKS=192
      ./xmlchange OCN_ROOTPE=1216

      ./xmlchange LND_NTASKS=1152
      ./xmlchange ROF_NTASKS=1152
      ./xmlchange ICE_NTASKS=256
      ./xmlchange LND_ROOTPE=256
      ./xmlchange ROF_ROOTPE=256

     elif [ "${nnodes}" == "22" ]; then

      ### Current defaults for M layout on chrysalis for F-Cases
      # The custom-22 was previously tested for bi-grid with L80-QBO

      ./xmlchange CPL_NTASKS=2816
      ./xmlchange ATM_NTASKS=2816   #let atm run on all CPU
      ./xmlchange OCN_NTASKS=384
      ./xmlchange OCN_ROOTPE=2432

      ./xmlchange LND_NTASKS=2304
      ./xmlchange ROF_NTASKS=2304
      ./xmlchange ICE_NTASKS=512
      ./xmlchange LND_ROOTPE=512
      ./xmlchange ROF_ROOTPE=512

     elif [ "${nnodes}" == "44" ]; then

      ### Current defaults for M layout on chrysalis for F-Cases
      # The custom-22 was previously tested for bi-grid with L80-QBO

      ./xmlchange CPL_NTASKS=5632
      ./xmlchange ATM_NTASKS=5632   #let atm run on all CPU
      ./xmlchange OCN_NTASKS=768
      ./xmlchange OCN_ROOTPE=4864

      ./xmlchange LND_NTASKS=4608
      ./xmlchange ROF_NTASKS=4608
      ./xmlchange ICE_NTASKS=1024
      ./xmlchange LND_ROOTPE=1024
      ./xmlchange ROF_ROOTPE=1024
 
     elif [ "${nnodes}" == "66" ]; then

      ### Current defaults for M layout on chrysalis for F-Cases
      # The custom-22 was previously tested for bi-grid with L80-QBO

      ./xmlchange CPL_NTASKS=8448
      ./xmlchange ATM_NTASKS=8448   #let atm run on all CPU
      ./xmlchange OCN_NTASKS=1152
      ./xmlchange OCN_ROOTPE=7296

      ./xmlchange LND_NTASKS=6912
      ./xmlchange ROF_NTASKS=6912
      ./xmlchange ICE_NTASKS=1536
      ./xmlchange LND_ROOTPE=1536
      ./xmlchange ROF_ROOTPE=1536

     elif [ "${nnodes}" == "80" ]; then
      	ncore=64	     
      ./xmlchange CPL_NTASKS=$(($nnodes * $ncore))
      ./xmlchange ATM_NTASKS=$(($nnodes * $ncore))
      ./xmlchange OCN_NTASKS=$(($nnodes * $ncore))
      ./xmlchange LND_NTASKS=$(($nnodes * $ncore))
      ./xmlchange ROF_NTASKS=$(($nnodes * $ncore))
      ./xmlchange ICE_NTASKS=$(($nnodes * $ncore))


     elif [ "${nnodes}" == "110" ]; then

      ### Current defaults for M layout on chrysalis for F-Cases
      # The custom-22 was previously tested for bi-grid with L80-QBO

      ./xmlchange CPL_NTASKS=14080
      ./xmlchange ATM_NTASKS=14080   #let atm run on all CPU
      ./xmlchange OCN_NTASKS=1920
      ./xmlchange OCN_ROOTPE=12160

      ./xmlchange LND_NTASKS=11520
      ./xmlchange ROF_NTASKS=11520
      ./xmlchange ICE_NTASKS=2560
      ./xmlchange LND_ROOTPE=2560
      ./xmlchange ROF_ROOTPE=2560
 
 
     else

       echo 'ERRROR: unsupported layout '${PELAYOUT}
       exit 401

    fi

    popd

fi

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
    local repo=E3SM

    echo "Cloning $repo repository branch $BRANCH under $path"
    if [ -d "${path}" ]; then
        echo "ERROR: Directory already exists. Not overwriting"
        exit 20
    fi
    mkdir -p ${path}
    pushd ${path}

    # This will put repository, with all code
    git clone git@github.com:E3SM-Project/${repo}.git .

    # Check out desired branch
    git checkout ${BRANCH}

    # Custom addition
    if [ "${CHERRY}" != "" ]; then
        echo ----- WARNING: adding git cherry-pick -----
        for commit in "${CHERRY[@]}"
        do
            echo ${commit}
            git cherry-pick ${commit}
        done
        echo -------------------------------------------
    fi

    # Bring in all submodule components
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

    if [[ ${PELAYOUT} == custom-* ]];
    then
        layout="M" # temporary placeholder for create_newcase
    else
        layout=${PELAYOUT}

    fi
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
        --pecount ${layout}

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

    # Setup some CIME directories
    ./xmlchange EXEROOT=${CASE_BUILD_DIR}
    ./xmlchange RUNDIR=${CASE_RUN_DIR}

    # Short term archiving
    ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING^^}
    ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

    # Build with COSP, except for a data atmosphere (datm)
    if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then 
      echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`


    # Custom user_nl
    user_nl

    ## Xue TEST8
    ./xmlchange CPL_SEQ_OPTION=RASM_OPTION2
    ./xmlchange LND2ATM_FMAPNAME='/lcrc/group/e3sm/data/inputdata/cpl/gridmaps/r025/map_r025_to_ne120pg2_traave.20240328.nc'
    ./xmlchange LND2ATM_SMAPNAME='/lcrc/group/e3sm/data/inputdata/cpl/gridmaps/r025/map_r025_to_ne120pg2_traave.20240328.nc'
    ./xmlchange LND_DOMAIN_FILE='domain.lnd.r025_RRSwISC6to18E3r5.240402.nc'
    ./xmlchange ROF2OCN_ICE_RMAPNAME='cpl/cpl6/map_r025_to_RRSwISC6to18E3r5.cstmnn.r250e1250_58NS.20240328.nc'
    # Finally, run CIME case.setup
    ./case.setup --reset
    # priority
    #./xmlchange --force JOB_QUEUE=priority

    popd
}

#-----------------------------------------------------
case_build() {

    pushd ${CASE_SCRIPTS_DIR}

    # do_case_build = false
    if [ "${do_case_build,,}" != "true" ]; then

        echo $'\n----- case_build -----\n'

        if [ "${OLD_EXECUTABLE}" == "" ]; then
            # Ues previously built executable, make sure it exists
            if [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
                echo 'Skipping build because $do_case_build = '${do_case_build}
            else
                echo 'ERROR: $do_case_build = '${do_case_build}' but no executable exists for this case.'
                exit 297
            fi
        else
            # If absolute pathname exists and is executable, reuse pre-exiting executable
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

    # do_case_build = true
    else

        echo $'\n----- Starting case_build -----\n'

        # Turn on debug compilation option if requested
        if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
            ./xmlchange DEBUG=${DEBUG_COMPILE^^}
        fi

        # Run CIME case.build
        ./case.build

    fi

    # Some user_nl settings won't be updated to *_in files under the run directory
    # Call preview_namelists to make sure *_in and user_nl files are consistent.
    echo $'\n----- Preview namelists -----\n'
    ./preview_namelists

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Set simulation start date
    ./xmlchange RUN_STARTDATE=${START_DATE}

    # Segment length
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}

    # Restart frequency
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}

    # Coupler history
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}

    # Coupler budgets (always on)
    ./xmlchange BUDGETS=TRUE

    # Set resubmissions
    if (( RESUBMIT > 0 )); then
        ./xmlchange RESUBMIT=${RESUBMIT}
    fi

    # Run type
    # Start from default of user-specified initial conditions
    if [ "${MODEL_START_TYPE,,}" == "initial" ]; then
        ./xmlchange RUN_TYPE="startup"
        ./xmlchange CONTINUE_RUN="FALSE"

    # Continue existing run
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

    # Patch mpas streams files
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
    
    # Run CIME case.submit
    ./case.submit

    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'

    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=`basename $0`
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp ${this_script_name} ${script_provenance_dir}/${script_provenance_name}

}

#-----------------------------------------------------
# Silent versions of popd and pushd
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

# Now, actually run the script
#-----------------------------------------------------
main

