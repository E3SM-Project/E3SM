#!/bin/bash -fe

# E3SM Water Cycle v2 run_e3sm script template.
#
# Inspired by v1 run_e3sm script as well as SCREAM group simplified run script.
#
# Bash coding style inspired by:
# http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

main() {

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=ruby
readonly PROJECT="cbronze"

# Simulation
readonly COMPSET="F20TR_chemUCI-Linozv3"
readonly RESOLUTION="CAx32v1pg2_CAx32v1pg2"
readonly CASE_NAME="20231003.t9.BRC006-CA.UCI6k.E3SM_RRMtrcdat.CArrm.amip.nudge.${MACHINE}"
# If this is part of a simulation campaign, ask your group lead about using a case_group label
# readonly CASE_GROUP=""

# Code and compilation
readonly CHECKOUT="20221110_CAfire"
readonly BRANCH="BRC-RRM" 
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="2020-09-05"

# Additional options for 'branch' and 'hybrid'
#readonly GET_REFCASE=TRUE
#readonly RUN_REFDIR="/global/cscratch1/sd/tang30/E3SM_simulations/tst.20221118.CArrm.amip.chemUCI_Linozv3.cori-knl/tests/custom-30_1x5_ndays/dt_remap_tstep4x_5days/archive/rest/2010-01-06-00000"
#readonly RUN_REFCASE="tst.20221118.CArrm.amip.chemUCI_Linozv3.cori-knl"
#readonly RUN_REFDATE="2010-01-06"   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

# Set paths
readonly CODE_ROOT="/p/lustre1/E3SMfire/ke2/E3SM_models/CAfire-Brc-RRMtrc/E3SM"
readonly CASE_ROOT="/p/lustre1/E3SMfire/ke2/E3SM_simulations/${CASE_NAME}"

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'S_2x5_ndays', 'M_1x10_ndays', 'M80_1x10_ndays'
#  or 'production' for full simulation
#readonly run='production'
readonly run='custom-20_ndays_UV50hU6hL-s20200905-trcdat-brc-decay'
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
  readonly WALLTIME="10:30:00"
  #readonly STOP_OPTION=${units}
  #readonly STOP_N=${length}
  #readonly REST_OPTION=${STOP_OPTION}
  #readonly REST_N=${STOP_N}
  readonly STOP_OPTION="ndays"
  readonly STOP_N="5"
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  #readonly RESUBMIT=${resubmit}
  readonly RESUBMIT=0
  readonly DO_SHORT_TERM_ARCHIVING=false

else

  # Production simulation
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="M"
  readonly WALLTIME="00:30:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="15"
  readonly REST_OPTION="nyears"
  readonly REST_N="1"
  readonly RESUBMIT="0"
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
 aer_drydep_list                = 'bc_a1', 'bc_a3', 'bc_a4', 'dst_a1', 'dst_a3', 'mom_a1', 'mom_a2', 'mom_a3', 'mom_a4', 'ncl_a1', 'ncl_a2',
         'ncl_a3', 'num_a1', 'num_a2', 'num_a3', 'num_a4', 'pom_a1', 'pom_a3', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3',
         'soa_a1', 'soa_a2', 'soa_a3', 'brc_a1', 'brc_a3', 'brc_a4'
 aer_wetdep_list                = 'bc_a1', 'bc_a3', 'bc_a4', 'dst_a1', 'dst_a3', 'mom_a1', 'mom_a2', 'mom_a3', 'mom_a4', 'ncl_a1', 'ncl_a2',
         'ncl_a3', 'num_a1', 'num_a2', 'num_a3', 'num_a4', 'pom_a1', 'pom_a3', 'pom_a4', 'so4_a1', 'so4_a2', 'so4_a3',
         'soa_a1', 'soa_a2', 'soa_a3', 'brc_a1', 'brc_a3', 'brc_a4'

 nhtfrq =   0,-6,-6,-24,-1
 mfilt  = 1,4,4,1,24
 avgflag_pertape = 'A','I','A','A','I'
 !fincl1 = 'CME','DCQ','DTCOND','DTCORE','TTGW','TTEND_TOT','PTTEND','PTEQ','PTECLDLIQ','PTECLDICE','WSUB','WLARGE','VOR','DIV','DYN_PS','DIV_Qflux'
 fincl1 = 'CME','DCQ','DTCOND','DTCORE','TTGW','brc_a4', 'brc_a3','brc_a1'
 fincl2 = 'U','V','T','Q','OMEGA500','OMEGA850','PS','U10','TMQ','Z500','TUQ','TVQ','Z3','PHIS','QREFHT','TREFHT','brc_a4', 'brc_a3','brc_a1'
 fincl3 = 'PRECT','TROP_P','TROP_Z'
 fincl4 = 'SAODVIS','AODVIS','Mass_bc','Mass_pom','MASS','PS','brc_a4', 'brc_a3','brc_a1','QRL','QRS','T','Z3','EXTINCT','OMEGA','TROP_P','TROP_Z','FSNT','FSNTC','FLNT','FLNTC','SWCF','LWCF','Mass_mom','Mass_ncl','Mass_soa','Mass_so4','Mass_dst'
 fincl5 = 'SAODVIS','AODVIS','Mass_bc','Mass_pom','MASS','PS','brc_a4', 'bc_a4','pom_a4','QRL','QRS','T','Z3','EXTINCT','OMEGA','TROP_P','TROP_Z','FSNT','FSNTC','FLNT','FLNTC','SWCF','LWCF','Mass_mom','Mass_ncl','Mass_soa','Mass_so4','Mass_dst','plume_height_EM','zmidr_ph','pmid_ph','tfld_ph','relhum_ph','qh2o_ph','ufld_ph','vfld_ph','SCO','TCO','O3','CO','bc_a4_XFRC','bc_a4_CLXF','brc_a4_XFRC','brc_a4_CLXF','pom_a4_XFRC','pom_a4_CLXF'
 !fincl5 = 'VOR:I','DIV:I','DYN_PS:I','DIV_Qflux:I'

 is_output_interactive_volc = .true.

 tropopause_e90_thrd		= 80.0e-9

 history_gaschmbudget_2D = .false.
 
 linoz_psc_t = 198.0
 rad_climate            = 'A:H2OLNZ:H2O', 'N:O2:O2', 'N:CO2:CO2',
         'A:O3:O3', 'A:N2OLNZ:N2O', 'A:CH4LNZ:CH4',
         'N:CFC11:CFC11', 'N:CFC12:CFC12', 
         'M:mam4_mode1:/p/lustre2/ke2/input_data/modes/modal_optics_mode1.nc',
          'M:mam4_mode2:/p/lustre2/ke2/input_data/modes/modal_optics_mode2.nc',
          'M:mam4_mode3:/p/lustre2/ke2/input_data/modes/modal_optics_mode3.nc',
          'M:mam4_mode4:/p/lustre2/ke2/input_data/modes/modal_optics_mode4.nc'

 sad_file		= '${input_data_dir}/atm/waccm/sulf/SAD_SULF_1849-2100_1.9x2.5_c090817.nc'
 
 !ncdata = '/p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/tst.20221118.CArrm.amip.chemUCI_Linozv3.cori-knl.eam.i.2010-02-01-00000.nc'
 !ncdata = '/p/lustre2/zhang73/E3SM_simulations/tst2A.20221213.CArrm.amip.nudge.quartz/tests/custom-30_ndays/run/tst2A.20221213.CArrm.amip.nudge.quartz.eam.i.2020-10-01-00000.nc'
 !ncdata = '/p/lustre1/E3SMfire/zhang73/E3SM_simulations/tst2B.20221213.CArrm.amip.nudge.quartz/tests/custom-30_ndays_UV5dU6hL/run/archive_r_i/tst2B.20221213.CArrm.amip.nudge.quartz.eam.i.2020-10-30-00000.nc'
!ncdata = '/p/lustre1/E3SMfire/zhang73/WRF_post/tign1_W2.0/UV50hU6hL-s20200901-E4.syrah.eam.i.2020-09-06-00000.wrf_tign1_W2.0_2020-09-06_04:30:00_CAx32v1np4_monotr.nc'
!ncdata = '/p/lustre1/E3SMfire/zhang73/WRF_post/tign1_W2.0/UV50hU6hL-s20200901-E4.syrah.eam.i.2020-09-06-00000.bc_brc_num_wrf.L72.tign1_W2.0_2020-09-06_04:30:00_CAx32v1np4_monotr.nc'
 !ncdata = '/p/lustre1/E3SMfire/zhang73/WRF_post/Creek4E3SM/UV50hU6hL-s20200901-E4.syrah.eam.i.2020-09-06-00000.bc_pom_num_wrf.L72.Creek4E3SM_2020-09-05_23:20:00_CAx32v1np4_monotr.nc'
ncdata='/p/lustre1/E3SMfire/zhang73/E3SM_simulations/tst2C.20221213.CArrmE4.amip.nudge.syrah/tests/custom-64_ndays_UVT50hU6hL-s20200901-E4/run/tst2C.20221213.CArrmE4.amip.nudge.syrah.eam.i.2020-09-05-00000.nc'
!ncdata='/p/lustre1/E3SMfire/zhang73/WRF_post/HRRR-Creek/UV50hU6hL-s20200901-E4.syrah.eam.i.2020-09-06-00000.nc'


 theta_hydrostatic_mode=.false.
 tstep_type		=  9
 deep_scheme = 'off'

 inithist = 'DAILY'

 Nudge_Model        =.true.
 Nudge_Path         ='/p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/L72.mono_fv2fv_CAx32v1pg2.UVTQ/'
 Nudge_File_Template='HICCUP.ERA5.%y-%m-%d-%s.72levs.mono_fv2fv_CAx32v1pg2.nc'
 Nudge_Times_Per_Day    =8   
 Nudge_Tau    		=50.
 Nudge_Tau2    		=6.
 Model_Times_Per_Day    =1152
 Nudge_Uprof            =2   
 Nudge_Ucoef            =0.5
 Nudge_Vprof            =2   
 Nudge_Vcoef            =0.5
 Nudge_Tprof            =0
 Nudge_Tcoef            =0.5
 Nudge_Qprof            =0
 Nudge_Qcoef            =0.5
 Nudge_PSprof           =0   
 Nudge_PScoef           =0.00
 Nudge_Beg_Year         =2019
 Nudge_Beg_Month        =10
 Nudge_Beg_Day          =1  
 Nudge_End_Year         =2021
 Nudge_End_Month        =12  
 Nudge_End_Day          =31  
 Nudge_Hwin_lo          =1.0 
 Nudge_Hwin_hi          =0.0 
 Nudge_Hwin_lat0        =37.2
 Nudge_Hwin_latWidth    =17.0
 Nudge_Hwin_latDelta    =1.4 
 Nudge_Hwin_lon0        =240.6
 Nudge_Hwin_lonWidth    =17.0
 Nudge_Hwin_lonDelta    =1.4
 Nudge_Vwin_lo          =0.0
 Nudge_Vwin_hi          =1.0
 Nudge_Vwin_Hindex      =100.0
 Nudge_Vwin_Hdelta      =10.
 Nudge_Vwin_Lindex      =0.0
 Nudge_Vwin_Ldelta      =0.1 !cannot.le.0
 Nudge_File_Ntime       =1
 Nudge_Method       ='Linear'
 Nudge_Loc_PhysOut  =.true.
 Nudge_CurrentStep  =.true.



!!Jinbo Xie added
chlorine_loading_file          = '/p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/Linoz_Chlorine_Loading_CMIP6_Hist_SSP370_0003-2503_c20210202.nc'
 ext_frc_specifier              = 'NO2         ->/p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_NO2_aircraft_vertical_1750-2025_1.9x2.5_c221208_sameafter2014.nc'
         'SO2         -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/elev_specifier/cmip6_ssp370_mam4_so2_elev_1850-2100_c221208.nc',
         'SOAG        -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/elev_specifier/cmip6_ssp370_mam4_soag_elev_1850-2100_c221208.nc',
         !'bc_a4       -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/elev_specifier/cmip6_ssp370_mam4_bc_a4_elev_1850-2100_c221208.nc',
         !'bc_a4       -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/CAx32v1pg2/cmip6_ssp370_mam4_bc_a4_elev_1850-2100_c221208.CAx32v1pg2.nc',
         'bc_a4       -> /p/lustre1/E3SMfire/ke2/wildfires/emis/cmip6_ssp370_UCI_6k_mam4_bc_a4_elev_2020_c221208.CAx32v1pg2_corrected.nc',
         !'bc_a4       -> /p/lustre1/E3SMfire/zhang73/WRF_post/tign1_W0.5/wrfout_d03-nt100_BB_bc_a4_elev_2020-09-06_CAx32v1pg2_highorder_c230417.nc',
         'num_a1      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/elev_specifier/cmip6_ssp370_mam4_num_a1_elev_1850-2100_c221208.nc',
         'num_a2      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/elev_specifier/cmip6_ssp370_mam4_num_a2_elev_1850-2100_c221208.nc',
         !'num_a4      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/elev_specifier/cmip6_ssp370_mam4_num_a4_elev_1850-2100_c221208.nc',
         'num_a4      -> /p/lustre1/E3SMfire/ke2/wildfires/emis/cmip6_ssp370_UCI_6k_mam4_num_a4_elev_2020_c221208.CAx32v1pg2_corrected.nc',
         'pom_a4      -> /p/lustre1/E3SMfire/ke2/wildfires/emis/cmip6_ssp370_mam4_pom_a4_elev_1850-2100_c221208.CAx32v1pg2.nc',
         'brc_a4      -> /p/lustre1/E3SMfire/ke2/wildfires/emis/cmip6_ssp370_UCI_6k_mam4_BRC-CA_a4_elev_2020_c221208.CAx32v1pg2_corrected.nc',
         'so4_a1      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/elev_specifier/cmip6_ssp370_mam4_so4_a1_elev_1850-2100_c221208.nc',
         'so4_a2      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/elev_specifier/cmip6_ssp370_mam4_so4_a2_elev_1850-2100_c221208.nc'
ext_frc_type           = 'SERIAL'
!ext_frc_cycle_yr       = 2014 
srf_emis_specifier             = 'C2H4      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_C2H4_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'C2H6      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_C2H6_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'C3H8      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_C3H8_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'CH2O      -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_CH2O_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'CH3CHO    -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_CH3CHO_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'CH3COCH3  -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_CH3COCH3_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'CO        -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_CO_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'DMS       -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc',
         'E90    -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions_E90_surface_1750-2025_1.9x2.5_c221208_sameafter2014.nc',
         'ISOP   -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_ISOP_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'NO     -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/emission/emissions-cmip6_e3sm_NO_surface_1850-2025_1.9x2.5_c221208_sameafter2014.nc',
         'SO2    -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/surf_specifier/cmip6_ssp370_mam4_so2_surf_1850-2100_c221208.nc',
         !'bc_a4  -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/surf_specifier/cmip6_ssp370_mam4_bc_a4_surf_1850-2100_c221208.nc',
         'bc_a4  -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/CAx32v1pg2/cmip6_ssp370_mam4_bc_a4_surf_1850-2100_c221208.CAx32v1pg2.nc',
         'num_a1 -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/surf_specifier/cmip6_ssp370_mam4_num_a1_surf_1850-2100_c221208.nc',
         'num_a2 -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/surf_specifier/cmip6_ssp370_mam4_num_a2_surf_1850-2100_c221208.nc',
         'num_a4 -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/surf_specifier/cmip6_ssp370_mam4_num_a4_surf_1850-2100_c221208.nc',
         'pom_a4 -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/surf_specifier/cmip6_ssp370_mam4_pom_a4_surf_1850-2100_c221208.nc',
         'so4_a1 -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/surf_specifier/cmip6_ssp370_mam4_so4_a1_surf_1850-2100_c221208.nc',
         'so4_a2 -> /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/surf_specifier/cmip6_ssp370_mam4_so4_a2_surf_1850-2100_c221208.nc',
 	 !
srf_emis_type          = 'SERIAL'
!srf_emis_cycle_yr      = 2014
	 linoz_data_file                = 'linv3_1849-2015_2010JPL_cmip6_historical_10deg_58km_c221208_sameafter2014.nc'
         linoz_data_path                = '/p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm'
	 !
&chem_surfvals_nl
  bndtvghg               = '/p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/atm/GHG_CMIP_SSP370-1-2-1_Annual_Global_0000-2500_c221208.nc'

 prescribed_volcaero_datapath           = ''
 prescribed_volcaero_file               = ''
 !prescribed_volcaero_filetype           = 'VOLC_CMIP6'
 !prescribed_volcaero_type               = 'CYCLICAL'
 !prescribed_volcaero_cycle_yr 		= 2015 

  mode_defs              = 'mam4_mode1:accum:=', 'A:num_a1:N:num_c1:num_mr:+',
         'A:so4_a1:N:so4_c1:sulfate:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:pom_a1:N:pom_c1:p-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:brc_a1:N:brc_c1:b-organic:/p/lustre2/ke2/input_data/modes/brcpho_rrtmg_c130709.nc:+',
         'A:soa_a1:N:soa_c1:s-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
         'A:bc_a1:N:bc_c1:black-c:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:dst_a1:N:dst_c1:dust:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+',
         'A:ncl_a1:N:ncl_c1:seasalt:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a1:N:mom_c1:m-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam4_mode2:aitken:=',
         'A:num_a2:N:num_c2:num_mr:+',
         'A:so4_a2:N:so4_c2:sulfate:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:soa_a2:N:soa_c2:s-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
         'A:ncl_a2:N:ncl_c2:seasalt:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a2:N:mom_c2:m-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam4_mode3:coarse:=',
         'A:num_a3:N:num_c3:num_mr:+',
         'A:dst_a3:N:dst_c3:dust:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+',
         'A:ncl_a3:N:ncl_c3:seasalt:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:so4_a3:N:so4_c3:sulfate:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:bc_a3:N:bc_c3:black-c:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:pom_a3:N:pom_c3:p-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:brc_a3:N:brc_c3:b-organic:/p/lustre2/ke2/input_data/modes/brcpho_rrtmg_c130709.nc:+',
         'A:soa_a3:N:soa_c3:s-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
         'A:mom_a3:N:mom_c3:m-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam4_mode4:primary_carbon:=', 'A:num_a4:N:num_c4:num_mr:+',
         'A:pom_a4:N:pom_c4:p-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:brc_a4:N:brc_c4:b-organic:/p/lustre2/ke2/input_data/modes/brcpho_rrtmg_c130709.nc:+',
         'A:bc_a4:N:bc_c4:black-c:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:mom_a4:N:mom_c4:m-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam4_mode1_nobrc:accum:=', 'A:num_a1:N:num_c1:num_mr:+',
         'A:so4_a1:N:so4_c1:sulfate:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:pom_a1:N:pom_c1:p-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:soa_a1:N:soa_c1:s-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
         'A:bc_a1:N:bc_c1:black-c:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:dst_a1:N:dst_c1:dust:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+',
         'A:ncl_a1:N:ncl_c1:seasalt:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a1:N:mom_c1:m-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam4_mode3_nobrc:coarse:=',
         'A:num_a3:N:num_c3:num_mr:+',
         'A:dst_a3:N:dst_c3:dust:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+',
         'A:ncl_a3:N:ncl_c3:seasalt:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:so4_a3:N:so4_c3:sulfate:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:bc_a3:N:bc_c3:black-c:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:pom_a3:N:pom_c3:p-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:soa_a3:N:soa_c3:s-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
         'A:mom_a3:N:mom_c3:m-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam4_mode4_nobrc:primary_carbon:=', 'A:num_a4:N:num_c4:num_mr:+',
         'A:pom_a4:N:pom_c4:p-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc:+',
         'A:bc_a4:N:bc_c4:black-c:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:mom_a4:N:mom_c4:m-organic:/usr/gdata/climdat/ccsm3data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc'

EOF

cat << EOF >> user_nl_elm
 check_finidat_year_consistency = .false.
 check_dynpft_consistency = .false.
 check_finidat_fsurdat_consistency = .false.
 !finidat = ' '
 !finidat = '/p/lustre1/E3SMfire/zhang73/E3SM_simulations/tst2B.20221213.CArrm.amip.nudge.quartz/tests/custom-30_ndays_UV5dU6hL/run/archive_r_i/tst2B.20221213.CArrm.amip.nudge.quartz.elm.r.2020-10-30-00000.nc'
 finidat = '/p/lustre1/E3SMfire/zhang73/E3SM_simulations/tst2C.20221213.CArrmE4.amip.nudge.syrah/tests/custom-64_ndays_UV50hU6hL-s20200901-E4/run/tst2C.20221213.CArrmE4.amip.nudge.syrah.elm.r.2020-09-07-00000.nc'
 !finidat = '/p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/tst.20221118.CArrm.amip.chemUCI_Linozv3.cori-knl.elm.r.2010-02-01-00000.nc'
 !fsurdat = '/global/cscratch1/sd/tang30/nudge/surfdata_CAx32v1pg2_simyr2020_c221211.nc'
EOF

cat <<EOF >> user_nl_cice
  model_year_align               = 1990
  stream_fldfilename             = '/p/lustre2/zhang73/forQi/E3SM_simulations/inputdata/sst_weekly_cdcunits_1x1_1990-01-01-2022-03-13_n3.nc'
  stream_fldvarname              = 'ice_cov'
  stream_year_first              = 1990
  stream_year_last               = 2022
EOF

cat <<EOF >> user_nl_docn
streams = "docn.streams.txt.prescribed 1990 1990 2022"
EOF

cat <<EOF >> user_docn.streams.txt.prescribed
<?xml version="1.0"?>
<file id="stream" version="1.0">
<dataSource>
   GENERIC
</dataSource>
<domainInfo>
  <variableNames>
     time    time
        xc      lon
        yc      lat
        area    area
        mask    mask
  </variableNames>
  <filePath>
     /usr/gdata/climdat/ccsm3data/inputdata/ocn/docn7
  </filePath>
  <fileNames>
     domain.ocn.1x1.111007.nc
  </fileNames>
</domainInfo>
<fieldInfo>
   <variableNames>
     SST_cpl t
   </variableNames>
   <filePath>
     /p/lustre2/zhang73/forQi/E3SM_simulations/inputdata
   </filePath>
   <fileNames>
     sst_weekly_cdcunits_1x1_1990-01-01-2022-03-13_n3.nc
   </fileNames>
   <offset>
      0
   </offset>
</fieldInfo>
</file>
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
    elif [ "${MACHINE}" == "cori-knl" ]; then
        #ncore=68
        # for CA RRM to use existing mpas partition file
        ncore=64
    elif [ "${MACHINE}" == "quartz" ]; then
	ncore=36
    elif [ "${MACHINE}" == "syrah" ]; then
	ncore=16
    elif [ "${MACHINE}" == "ruby" ]; then
	ncore=56
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
    local usr_mech_infile="${CODE_ROOT}/components/eam/chem_proc/inputs/pp_chemUCI_linozv3_mam4_brc_decay_resus_mom_soag_tag.in"
    #echo '[QT] Changing chemistry to :'${usr_mech_infile}
    #./xmlchange --id CAM_CONFIG_OPTS --append --val='-usr_mech_infile '${usr_mech_infile}
    #./xmlchange --id CAM_CONFIG_OPTS --append --val="-cosp -chem superfast_mam4_brc_resus_mom  -usr_mech_infile ${usr_mech_infile}"
    ./xmlchange --id CAM_CONFIG_OPTS --append --val="-chem superfast_mam4_brc_resus_mom  -usr_mech_infile ${usr_mech_infile}"
    # This is specific to the CAx32 RRM grid
    ./xmlchange EPS_AGRID=1e-9

    #./xmlchange PIO_VERSION=1

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

