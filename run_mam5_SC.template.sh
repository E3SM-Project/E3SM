#!/bin/bash -fe

# E3SM Water Cycle v2 run_e3sm script template.
#
# Inspired by v1 run_e3sm script as well as SCREAM group simplified run script.
#
# Bash coding style inspired by:
# http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

# TO DO:
# - custom pelayout

main() {

# For debugging, uncomment line below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=cori-knl
readonly PROJECT="e3sm"

# Simulation
readonly COMPSET="F20TR"
readonly RESOLUTION="ne30pg2_EC30to60E2r2"
# BEFORE RUNNING : CHANGE the following CASE_NAME to desired value
readonly CASE_NAME="run.mam5_SC.maint-2.0.test"
# If this is part of a simulation campaign, ask your group lead about using a case_group label
# readonly CASE_GROUP=""

# Code and compilation
readonly CHECKOUT="E3SMv2-mam5/E3SM"
readonly BRANCH="maint-2.0-mam5" # this has no impact on setup, notation purpose only
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="0001-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
readonly RUN_REFDIR="/global/cscratch1/sd/forsyth/E3SMv2/v2.LR.piControl/init"
readonly RUN_REFCASE="20210625.v2rc3c-GWD.piControl.ne30pg2_EC30to60E2r2.chrysalis"
readonly RUN_REFDATE="1001-01-01"   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

# Set paths
# /global/homes/z/zke/E3SM_models/E3SMv2-mam5/E3SM
readonly CODE_ROOT="${HOME}/E3SM_models/${CHECKOUT}"
readonly CASE_ROOT="/global/cscratch1/sd/${USER}/E3SM_simulations/E3SMv2/${CASE_NAME}"

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays',
#               'M_1x10_ndays', 'M2_1x10_ndays', 'M80_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation
readonly run='XS_2x5_ndays'
if [ "${run}" != "production" ]; then

  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}

  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/tests/${run}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/tests/${run}/run
  readonly PELAYOUT=${layout}
  readonly WALLTIME="2:00:00"
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false

else

  # Production simulation
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="L"
  readonly WALLTIME="34:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="50"
  readonly REST_OPTION="nyears"
  readonly REST_N="5"
  readonly RESUBMIT="9"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history
readonly HIST_OPTION="nyears"
readonly HIST_N="5"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

# --- Now, do the work ---

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

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
 nhtfrq =   0,-24,-6,-6,-3,-24,0
 mfilt  = 1,30,120,120,240,30,1
 avgflag_pertape = 'A','A','I','A','A','A','I'
 fexcl1 = 'CFAD_SR532_CAL', 'LINOZ_DO3', 'LINOZ_DO3_PSC', 'LINOZ_O3CLIM', 'LINOZ_O3COL', 'LINOZ_SSO3', 'hstobie_linoz'
 fincl1 = 'extinct_sw_inp','extinct_lw_bnd7','extinct_lw_inp','CLD_CAL', 'TREFMNAV', 'TREFMXAV'
 fincl2 = 'FLUT','PRECT','U200','V200','U850','V850','Z500','OMEGA500','UBOT','VBOT','TREFHT','TREFHTMN:M','TREFHTMX:X','QREFHT','TS','PS','TMQ','TUQ','TVQ','TOZ', 'FLDS', 'FLNS', 'FSDS', 'FSNS', 'SHFLX', 'LHFLX', 'TGCLDCWP', 'TGCLDIWP', 'TGCLDLWP', 'CLDTOT', 'T250', 'T200', 'T150', 'T100', 'T050', 'T025', 'T010', 'T005', 'T002', 'T001', 'TTOP', 'U250', 'U150', 'U100', 'U050', 'U025', 'U010', 'U005', 'U002', 'U001', 'UTOP', 'FSNT', 'FLNT'
 fincl3 = 'PSL','T200','T500','U850','V850','UBOT','VBOT','TREFHT', 'Z700', 'TBOT:M'
 fincl4 = 'FLUT','U200','U850','PRECT','OMEGA500'
 fincl5 = 'PRECT','PRECC','TUQ','TVQ','QFLX','SHFLX','U90M','V90M'
 fincl6 = 'CLDTOT_ISCCP','MEANCLDALB_ISCCP','MEANTAU_ISCCP','MEANPTOP_ISCCP','MEANTB_ISCCP','CLDTOT_CAL','CLDTOT_CAL_LIQ','CLDTOT_CAL_ICE','CLDTOT_CAL_UN','CLDHGH_CAL','CLDHGH_CAL_LIQ','CLDHGH_CAL_ICE','CLDHGH_CAL_UN','CLDMED_CAL','CLDMED_CAL_LIQ','CLDMED_CAL_ICE','CLDMED_CAL_UN','CLDLOW_CAL','CLDLOW_CAL_LIQ','CLDLOW_CAL_ICE','CLDLOW_CAL_UN'
 fincl7 = 'O3', 'PS', 'TROP_P'

! Additional retuning
 clubb_tk1 = 268.15D0
 gw_convect_hcf = 10.0

!MAM5 added
&aerosol_nl
 seasalt_emis_scale             =  0.6
 sol_factb_interstitial         = 0.1D0
 sol_facti_cloud_borne          = 1.0D0
 sol_factic_interstitial                = 0.4D0
 sscav_tuning           = .true.
/


docosp    = .true.,
cosp_lite = .true.,
cosp_ncolumns        = 10
cosp_nradsteps       = 3
cosp_lmisr_sim       = .true.
cosp_lisccp_sim      = .true.
cosp_lmodis_sim      = .true.
cosp_llidar_sim      = .true.
history_amwg         = .true.
history_aero_optics  = .true.
history_aerosol      = .true.
history_budget       = .true.
history_verbose      = .true.
do_aerocom_ind3      = .true.


   mode_defs              = 'mam5_mode1:accum:=', 'A:num_a1:N:num_c1:num_mr:+',
         'A:so4_a1:N:so4_c1:sulfate:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+', 'A:pom_a1:N:pom_c1:p-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:soa_a1:N:soa_c1:s-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
 'A:bc_a1:N:bc_c1:black-c:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:dst_a1:N:dst_c1:dust:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+', 'A:ncl_a1:N:ncl_c1:seasalt:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a1:N:mom_c1:m-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
'mam5_mode2:aitken:=',
         'A:num_a2:N:num_c2:num_mr:+', 'A:so4_a2:N:so4_c2:sulfate:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:soa_a2:N:soa_c2:s-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
'A:ncl_a2:N:ncl_c2:seasalt:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a2:N:mom_c2:m-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
'mam5_mode3:coarse:=',
         'A:num_a3:N:num_c3:num_mr:+', 'A:dst_a3:N:dst_c3:dust:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+',
         'A:ncl_a3:N:ncl_c3:seasalt:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
 'A:so4_a3:N:so4_c3:sulfate:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:bc_a3:N:bc_c3:black-c:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
 'A:pom_a3:N:pom_c3:p-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:soa_a3:N:soa_c3:s-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 'A:mom_a3:N:mom_c3:m-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
 'mam5_mode4:primary_carbon:=', 'A:num_a4:N:num_c4:num_mr:+',
         'A:pom_a4:N:pom_c4:p-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
'A:bc_a4:N:bc_c4:black-c:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:mom_a4:N:mom_c4:m-organic:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
  'mam5_mode5:strat_coarse:=', 'A:num_a5:N:num_c5:num_mr:+',
         'A:so4_a5:N:so4_c5:sulfate:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc',

 rad_climate            = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2',
         'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4',
         'N:CFC11:CFC11', 'N:CFC12:CFC12',
         'M:mam5_mode1:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode1_rrtmg_c130628.nc',
         'M:mam5_mode2:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode2_rrtmg_c130628.nc',
         'M:mam5_mode3:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode3_rrtmg_c130628.nc',
         'M:mam5_mode4:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode4_rrtmg_c130628.nc',
         'M:mam5_mode5:/global/cscratch1/sd/zke/CESM_input/modes/mam4_mode3_rrtmg_aeronetdust_sig1.2_dgnl.40_c150219_ke07.nc',
 rad_diag_1 =  'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2',
         'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4',
         'N:CFC11:CFC11', 'N:CFC12:CFC12',

 rad_diag_2 =  'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2',
         'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4',
         'N:CFC11:CFC11', 'N:CFC12:CFC12',
          'M:mam5_mode1:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode1_rrtmg_c130628.nc',
         'M:mam5_mode2:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode2_rrtmg_c130628.nc',
         'M:mam5_mode3:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode3_rrtmg_c130628.nc',
         'M:mam5_mode4:/project/projectdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode4_rrtmg_c130628.nc'

&dust_nl
 dust_emis_fact         =  1.50D0
 soil_erod_file         = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/dst/dst_1.9x2.5_c090203.nc'
/

&chem_inparm
 chlorine_loading_file          = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart/ub/Linoz_Chlorine_Loading_CMIP6_0003-2017_c20171114.nc'
 chlorine_loading_type          = 'SERIAL'
 clim_soilw_file                = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart/dvel/clim_soilw.nc'
 depvel_file            = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart/dvel/depvel_monthly.nc'
 depvel_lnd_file                = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart/dvel/regrid_vegetation.nc'
 drydep_srf_file                = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mam/atmsrf_ne30pg2_200129.nc'
 exo_coldens_file               = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart/phot/exo_coldens.nc'
 ext_frc_specifier              =
         'SOAG        -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_soag_elev_1850-2014_c180205.nc',
         'SO2         -> /global/project/projectdirs/e3sm/zke/CESM_input/emis/vol_emis/cmip6_mam4_so2_elev_1850-2014_c180205_kzm_1990s_volcano_adjusted.nc',
         'bc_a4       -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_elev_1850-2014_c180205.nc',
         'num_a1      -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_elev_1850-2014_c180205.nc',
         'num_a2      -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_elev_1850-2014_c180205.nc',
         'num_a4      -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_elev_1850-2014_c180205.nc',
         'pom_a4      -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_elev_1850-2014_c180205.nc',
         'so4_a1      -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_elev_1850-2014_c180205.nc',
         'so4_a2      -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_elev_1850-2014_c180205.nc'
 ext_frc_type           = 'INTERP_MISSING_MONTHS'
 fstrat_list            = ' '
 linoz_data_file                = 'linoz1850-2015_2010JPL_CMIP6_10deg_58km_c20171109.nc'
 linoz_data_path                = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart/ub'
 linoz_data_type                = 'INTERP_MISSING_MONTHS'
 rsf_file               = '/global/cfs/cdirs/e3sm/inputdata/atm/waccm/phot/RSF_GT200nm_v3.0_c080811.nc'
 sad_file               = '/global/cfs/cdirs/e3sm/inputdata/atm/waccm/sulf/SAD_SULF_1950-2011_1.9x2.5_c130102.nc'
 sad_type               = 'SERIAL'
 season_wes_file                = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart/dvel/season_wes.nc'
 srf_emis_specifier             = 'DMS       -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc',
         'SO2       -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_surf_1850-2014_c180205.nc',
         'bc_a4     -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_surf_1850-2014_c180205.nc',
         'num_a1    -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_surf_1850-2014_c180205.nc',
         'num_a2    -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_surf_1850-2014_c180205.nc',
         'num_a4    -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_surf_1850-2014_c180205.nc',
         'pom_a4    -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_surf_1850-2014_c180205.nc',
         'so4_a1    -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_surf_1850-2014_c180205.nc',
         'so4_a2    -> /project/projectdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_surf_1850-2014_c180205.nc'
 srf_emis_type          = 'INTERP_MISSING_MONTHS'
 tracer_cnst_datapath           = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/oxid'
 tracer_cnst_file               = 'oxid_1.9x2.5_L26_1850-2015_c20181106.nc'
 tracer_cnst_filelist           = ''
 tracer_cnst_specifier          = 'cnst_O3:O3','OH','NO3','HO2'
 tracer_cnst_type               = 'INTERP_MISSING_MONTHS'
  xactive_prates         = .false.
 xs_coef_file           = '/global/cfs/cdirs/e3sm/inputdata/atm/waccm/phot/effxstex.txt'
 xs_long_file           = '/global/cfs/cdirs/e3sm/inputdata/atm/waccm/phot/temp_prs_GT200nm_JPL10_c130206.nc'
 xs_short_file          = '/global/cfs/cdirs/e3sm/inputdata/atm/waccm/phot/xs_short_jpl10_c130206.nc'
/

&chem_surfvals_nl
 bndtvghg               = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/ggas/GHG_CMIP-1-2-0_Annual_Global_0000-2014_c20180105.nc'
 co2vmr         = 0.000001e-6
 flbc_list              = ' '
 scenario_ghg           = 'RAMPED'
/


EOF

cat << EOF >> user_nl_elm
 hist_dov2xy = .true.,.true.
 hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA'
 hist_mfilt = 1,365
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'
EOF

cat << EOF >> user_nl_mosart
 rtmhist_fincl2 = 'RIVER_DISCHARGE_OVER_LAND_LIQ'
 rtmhist_mfilt = 1,365
 rtmhist_ndens = 2
 rtmhist_nhtfrq = 0,-24
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

    # This will put repository, with all code
    git clone git@github.com:E3SM-Project/${repo}.git .

    # Setup git hooks
    rm -rf .git/hooks
    git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
    git config commit.template .git/hooks/commit.template

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
			--pecount ${PELAYOUT}
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
			--pecount ${PELAYOUT}
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
    ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp -chem linoz_mam5_resus_mom_soag  -usr_mech_infile ' 
    ./xmlchange --id CAM_CONFIG_OPTS --append --val=${CODE_ROOT}'/components/eam/chem_proc/inputs/chem_mech_MAM5_SC.in'
    # Custom user_nl
    user_nl

    # Finally, run CIME case.setup
    ./case.setup --reset

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

        # Some user_nl settings won't be updated to *_in files under the run directory
        # Call preview_namelists to make sure *_in and user_nl files are consistent.
        ./preview_namelists

    fi

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
