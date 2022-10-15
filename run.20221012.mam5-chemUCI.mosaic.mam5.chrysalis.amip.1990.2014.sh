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

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=chrysalis
readonly PROJECT="e3sm"

# Simulation
readonly COMPSET="F20TR_chemUCI-Linozv3"
readonly RESOLUTION="ne30pg2_EC30to60E2r2"
readonly CASE_NAME="20221007.M5.AMIP.NGD-v3atm.M5-UCI-MOSAIC.1985.2014"
#readonly CASE_GROUP="v2.LR.chemUCI"

# Code and compilation
readonly CHECKOUT="1012-NGD-v3atm"
#readonly BRANCH="master"   #420e251265928a4122a5eb5d4eab4ae8110f84f2 05/172022  
readonly BRANCH="mam5_BFB_fix"     
#readonly CHERRY=( "7e4d1c9fec40ce1cf2c272d671f5d9111fa4dea7" "a5b1d42d7cd24924d0dbda95e24ad8d4556d93f1" ) # PR4349
readonly DEBUG_COMPILE=false

# Run options
readonly MODEL_START_TYPE="hybrid"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="1985-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
readonly RUN_REFDIR="/lcrc/group/e3sm/ac.zke/E3SM_inputs/rest/v2.LR.amip_0101/rest/1985-01-01-00000/"
readonly RUN_REFCASE="v2.LR.amip_0101"
readonly RUN_REFDATE="1985-01-01"   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

# Set paths
readonly CODE_ROOT="${HOME}/E3SM_models/${CHECKOUT}/v3atm"
#readonly CODE_ROOT="/qfs/people/tang338/E3SM_code/${CHECKOUT}"
readonly CASE_ROOT="/lcrc/group/e3sm/${USER}/E3SM_simulations/${CASE_NAME}"

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays', 
#               'M_1x10_ndays', 'M2_1x10_ndays', 'M80_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation
readonly run='production'
#readonly run='S_1x5_ndays'
if [ "${run}" != "production" ]; then

  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}

  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/tests/${run}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/tests/${run}/run
  #readonly PELAYOUT=${layout}
  readonly PELAYOUT="custom-10"
  readonly WALLTIME="00:20:00"
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
  #readonly PELAYOUT="M"
  readonly PELAYOUT="custom-30"
  readonly WALLTIME="15:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="10"
  readonly REST_OPTION="nyears"
  readonly REST_N="1"
  readonly RESUBMIT="2"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler histfalse 
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

readonly new_emis_dir="/compyfs/wumi635/inputdata/cam/chem/emis/CMIP6_emissions_1750_2015_2deg_FINAL"

cat << EOF >> user_nl_eam


&aerosol_nl
 seasalt_emis_scale             =  0.6
 sol_factb_interstitial         = 0.1D0
 sol_facti_cloud_borne          = 1.0D0
 sol_factic_interstitial                = 0.4D0
 sscav_tuning           = .true.
/

nhtfrq = 0, -24, -6,-6,-3,-24,0
 mfilt  = 1, 30,120,120,240,30,1
 avgflag_pertape = 'A','A','I','A','A','A','I'
 fexcl1 = 'CFAD_SR532_CAL', 'LINOZ_DO3', 'LINOZ_DO3_PSC', 'LINOZ_O3CLIM', 'LINOZ_O3COL', 'LINOZ_SSO3', 'hstobie_linoz'
 fincl1 = 'PS', 'angstrm', 'cod', 'cdr', 'cdnc', 'cdnum', 'icnum', 'clt', 'lcc', 'lwp', 'iwp', 'icc',
         'icnc', 'icr', 'LHFLX', 'SHFLX', 'OMEGA500', 'rh700', 'colrv', 'ccn', 'ccn.1bl', 'ccn.3bl', 'ptop', 'ttop',
         'rwp', 'lwp2', 'iwp2', 'autoconv', 'accretn', 'FSUTOA_d1', 'FSUTOAC_d1', 'FSUTOA', 'FSUTOAC', 'FLUTC', 'FLUT', 'PRECC',
         'PRECL', 'PRECT', 'TH7001000',
           'so4_a1','so4_a2','so4_a3','so4_a5','so4_c5','PS','num_a1','num_a2','num_a3','num_a4','num_a5',
          'SO2','Z3','AREA','bc_a1','bc_a4','dst_a1','dst_a3','pom_a1','pom_a4','soa_a1','soa_a2','soa_a3',
          'SO2','SO2_CLXF', 'SO2_XFRC',
 'dgnd_a01', 'dgnd_a02','dgnd_a03','dgnd_a04','dgnd_a05','dgnumwet1','dgnumwet2','dgnumwet3','dgnumwet4','dgnumwet5',
 'AODMODE1','AODMODE2','AODMODE3','AODMODE4','AODMODE5','NIMIX_IMM','NIMIX_CNT','NIMIX_DEP','NIHF', 'NIIMM','NIDEP','CCN3','extinct_sw_inp','extinct_lw_bnd7','extinct_lw_inp','CLD_CAL', 'TREFMNAV', 'TREFMXAV'

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

prescribed_volcaero_file = ‘’
prescribed_volcaero_datapath = ‘’

&chem_inparm
 chlorine_loading_file          = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/ub/Linoz_Chlorine_Loading_CMIP6_0003-2017_c20171114.nc'
 chlorine_loading_type          = 'SERIAL'
 clim_soilw_file                = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/dvel/clim_soilw.nc'
 depvel_file            = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/dvel/depvel_monthly.nc'
 depvel_lnd_file                = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/dvel/regrid_vegetation.nc'
 drydep_srf_file                = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mam/atmsrf_ne30pg2_200129.nc'
 exo_coldens_file               = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/phot/exo_coldens.nc'
 ext_frc_specifier              = 'NO2         -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_NO2_aircraft_vertical_1750-2015_1.9x2.5_c20170608.nc',
         'SO2         -> /lcrc/group/e3sm/ac.zke/E3SM_inputs/vol_emis/cmip6_mam4_so2_elev_1850-2014_c180205_kzm_1990s_volcano_adjusted.nc',
         'SOAG        -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_soag_elev_1850-2014_c180205.nc',
         'bc_a4       -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_elev_1850-2014_c180205.nc',
         'num_a1      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_elev_1850-2014_c180205.nc',
         'num_a2      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_elev_1850-2014_c180205.nc',
         'num_a4      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_elev_1850-2014_c180205.nc',
         'pom_a4      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_elev_1850-2014_c180205.nc',
         'so4_a1      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_elev_1850-2014_c180205.nc',
         'so4_a2      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_elev_1850-2014_c180205.nc'
 ext_frc_type           = 'INTERP_MISSING_MONTHS'
 fstrat_file            = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/ub/ubvals_b40.20th.track1_1996-2005_c110315.nc'
 fstrat_list            = ' '
 linoz_data_file                = 'linv3_1849-2015_2010JPL_cmip6_historical_10deg_58km_c20210625.nc'
 linoz_data_path                = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/ub'
 linoz_data_type                = 'INTERP_MISSING_MONTHS'
 o2_xsect_file          = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/phot/o2src.nc'
 rsf_file               = '/lcrc/group/e3sm/data/inputdata/atm/waccm/phot/RSF_GT200nm_v3.0_c080811.nc'
 sad_file               = '/lcrc/group/e3sm/data/inputdata/atm/waccm/sulf/SAD_SULF_1950-2011_1.9x2.5_c130102.nc'
 sad_type               = 'SERIAL'
 season_wes_file                = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/dvel/season_wes.nc'
 srf_emis_specifier             =
         'C2H4      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C2H4_surface_1850-2014_1.9x2.5_c20210323.nc',
         'C2H6      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C2H6_surface_1850-2014_1.9x2.5_c20210323.nc',
         'C3H8      -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C3H8_surface_1850-2014_1.9x2.5_c20210323.nc',
         'CH2O   -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH2O_surface_1850-2014_1.9x2.5_c20210323.nc',
         'CH3CHO    -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH3CHO_surface_1850-2014_1.9x2.5_c20210323.nc',
         'CH3COCH3  -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH3COCH3_surface_1850-2014_1.9x2.5_c20210323.nc',
         'CO     -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CO_surface_1850-2014_1.9x2.5_c20210323.nc',
         'DMS       -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc',
         'E90       -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions_E90_surface_1750-2015_1.9x2.5_c20210408.nc',
         'ISOP   -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_ISOP_surface_1850-2014_1.9x2.5_c20210323.nc',
         'NO     -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_NO_surface_1850-2014_1.9x2.5_c20210323.nc',
         'SO2       -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_surf_1850-2014_c180205.nc',
         'bc_a4     -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_surf_1850-2014_c180205.nc',
         'num_a1    -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_surf_1850-2014_c180205.nc',
         'num_a2    -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_surf_1850-2014_c180205.nc',
         'num_a4    -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_surf_1850-2014_c180205.nc',
         'pom_a4    -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_surf_1850-2014_c180205.nc',
         'so4_a1    -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_surf_1850-2014_c180205.nc',
         'so4_a2    -> /lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_surf_1850-2014_c180205.nc'
 srf_emis_type          = 'INTERP_MISSING_MONTHS'
 sulf_file              = '/lcrc/group/e3sm/data/inputdata/atm/waccm/sulf/sulfate.ar5_camchem_c130304.nc'
 tracer_cnst_cycle_yr           = 1995
 tracer_cnst_datapath           = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/methane/'
 tracer_cnst_file               = 'ch4_oxid_1.9x2.5_L26_1990-1999clim.c090804.nc'
 tracer_cnst_filelist           = ''
 tracer_cnst_specifier          = 'CH4','cnst_NO3:NO3', 'cnst_OH:OH'
 tracer_cnst_type               = 'CYCLICAL'
 xactive_prates         = .false.
 xs_long_file           = '/lcrc/group/e3sm/data/inputdata/atm/waccm/phot/temp_prs_GT200nm_JPL10_c130206.nc'
/

&chem_surfvals_nl
 bndtvghg               = '/lcrc/group/e3sm/data/inputdata/atm/cam/ggas/GHG_CMIP-1-2-0_Annual_Global_0000-2014_c20180105.nc'
 co2vmr         = 0.000001e-6
 scenario_ghg           = 'RAMPED'
/

&dust_nl
 dust_emis_fact         =  1.50D0
 soil_erod_file         = '/lcrc/group/e3sm/data/inputdata/atm/cam/dst/dst_1.9x2.5_c090203.nc'
/


&rad_cnst_nl
 icecldoptics           = 'mitchell'
 iceopticsfile          = '/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/iceoptics_c080917.nc'
 liqcldoptics           = 'gammadist'
 liqopticsfile          = '/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/F_nwvl200_mu20_lam50_res64_t298_c080428.nc'
mode_defs              =
         'mam5_mode1:accum:=', 'A:num_a1:N:num_c1:num_mr:+',
         'A:so4_a1:N:so4_c1:sulfate:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+', 'A:pom_a1:N:pom_c1:p-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:soa_a1:N:soa_c1:s-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 'A:bc_a1:N:bc_c1:black-c:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:dst_a1:N:dst_c1:dust:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/dust_aeronet_rrtmg_c141106.nc:+', 'A:ncl_a1:N:ncl_c1:seasalt:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a1:N:mom_c1:m-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam5_mode2:aitken:=',
         'A:num_a2:N:num_c2:num_mr:+', 'A:so4_a2:N:so4_c2:sulfate:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:soa_a2:N:soa_c2:s-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 'A:ncl_a2:N:ncl_c2:seasalt:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a2:N:mom_c2:m-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam5_mode3:coarse:=',
         'A:num_a3:N:num_c3:num_mr:+', 'A:dst_a3:N:dst_c3:dust:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/dust_aeronet_rrtmg_c141106.nc:+',
         'A:ncl_a3:N:ncl_c3:seasalt:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+', 'A:so4_a3:N:so4_c3:sulfate:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:bc_a3:N:bc_c3:black-c:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+', 'A:pom_a3:N:pom_c3:p-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:soa_a3:N:soa_c3:s-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 'A:mom_a3:N:mom_c3:m-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam5_mode4:primary_carbon:=', 'A:num_a4:N:num_c4:num_mr:+',
         'A:pom_a4:N:pom_c4:p-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+', 'A:bc_a4:N:bc_c4:black-c:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:mom_a4:N:mom_c4:m-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam5_mode5:strat_coarse:=', 'A:num_a5:N:num_c5:num_mr:+',
         'A:so4_a5:N:so4_c5:sulfate:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc',
 rad_climate            = 'A:H2OLNZ:H2O', 'N:O2:O2', 'N:CO2:CO2',
         'A:O3:O3', 'A:N2OLNZ:N2O', 'A:CH4LNZ:CH4',
         'N:CFC11:CFC11', 'N:CFC12:CFC12',
         'M:mam5_mode1:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode1_rrtmg_aeronetdust_c141106.nc',
         'M:mam5_mode2:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode2_rrtmg_c130628.nc',
         'M:mam5_mode3:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode3_rrtmg_aeronetdust_c141106.nc',
         'M:mam5_mode4:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode4_rrtmg_c130628.nc',
         'M:mam5_mode5:/lcrc/group/e3sm/ac.zke/E3SM_inputs/physprops/mam4_mode3_rrtmg_aeronetdust_sig1.2_dgnl.40_c150219_ke.nc'
 rad_diag_1 = 'A:H2OLNZ:H2O', 'N:O2:O2', 'N:CO2:CO2',
         'A:O3:O3', 'A:N2OLNZ:N2O', 'A:CH4LNZ:CH4',
         'N:CFC11:CFC11', 'N:CFC12:CFC12',
/



&tropopause_nl
 tropopause_climo_file          = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/ub/clim_p_trop.nc'
 tropopause_e90_thrd            = 80.0e-9
 tropopause_output_all          = .true.
/
 !history_gaschmbudget = .true.
 history_gaschmbudget_2D = .true.
 history_gaschmbudget_2D_levels = .true.
 history_UCIgaschmbudget_2D = .true.
 history_UCIgaschmbudget_2D_levels = .true.
 gaschmbudget_2D_L1_s =  1
 gaschmbudget_2D_L1_e = 26
 gaschmbudget_2D_L2_s = 27
 gaschmbudget_2D_L2_e = 38
 gaschmbudget_2D_L3_s = 39
 gaschmbudget_2D_L3_e = 58
 gaschmbudget_2D_L4_s = 59
 gaschmbudget_2D_L4_e = 72

 linoz_psc_t = 198.0

EOF

cat << EOF >> user_nl_elm
 check_finidat_year_consistency = .false.
 check_dynpft_consistency = .false.
 fsurdat = "${input_data_dir}/lnd/clm2/surfdata_map/surfdata_ne30pg2_simyr1850_c210402.nc"

EOF


}

# =====================================
# Customize MPAS stream files if needed
# =====================================


patch_mpas_streams() {

echo

}




# =====================================================
# Custom PE layout: custom-N where N is number of nodes
# =====================================================

custom_pelayout() {

if [[ ${PELAYOUT} == custom-* ]];
then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    # Number of cores per node (machine specific)
    if [ "${MACHINE}" == "chrysalis" ]; then
        ncore=64
    elif [ "${MACHINE}" == "compy" ]; then
        ncore=40
    else
        echo 'ERROR: MACHINE = '${MACHINE}' is not supported for custom PE layout.' 
        exit 400
    fi

    # Extract number of nodes
    tmp=($(echo ${PELAYOUT} | tr "-" " "))
    nnodes=${tmp[1]}

    # Customize
    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange NTASKS=$(( $nnodes * $ncore ))
    ./xmlchange NTHRDS=1
    ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
    ./xmlchange MAX_TASKS_PER_NODE=$ncore
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

    # Bring in all submodule components
    git submodule update --init --recursive

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

    # QT turn off cosp for testing
    ## Build with COSP, except for a data atmosphere (datm)
    #if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then 
    #  echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    #else
    #  echo $'\nConfiguring E3SM to use the COSP simulator\n'
    #  ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    #fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # QT changing chemistry mechanism
    #local usr_mech_infile="${CODE_ROOT}/components/eam/chem_proc/inputs/pp_chemUCI_linozv3_mam4_resus_mom_soag_tag.in"
    #echo '[QT] Changing chemistry to :'${usr_mech_infile}
    #./xmlchange --id CAM_CONFIG_OPTS --append --val='-usr_mech_infile '${usr_mech_infile}
    if [ "${run}" != "production" ]; then
       local CASE_ROOT_1="/lcrc/group/e3sm/${USER}/E3SM_simulations/${CASE_NAME}/tests/S_1x5_ndays/case_scripts/"      
    else
       local CASE_ROOT_1="$CASE_ROOT/case_scripts/"	    
    fi
    echo $CASE_ROOT_1
    local usr_mech_infile="${CODE_ROOT}/components/eam/chem_proc/inputs/pp_chemUCI_linozv3_mam5_resus_mom_soag_tag_cnst_no3_oh.in" 
    ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp -chem linoz_mam5_resus_mom_soag  -usr_mech_infile '${usr_mech_infile}
    # Custom user_nl
    user_nl
    ./xmlchange --id SSTICE_DATA_FILENAME --val='$DIN_LOC_ROOT/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc'

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

