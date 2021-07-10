#!/bin/bash

# E3SM Water Cycle v2 run_e3sm script template.
#
# Inspired by v1 run_e3sm script as well as SCREAM group simplified run script.
#
# Bash coding style inspired by:
# http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

# TO DO:
# - branch, hybrid restart
# - custom pelayout

main() {

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=cori-knl
readonly PROJECT="e3sm"

# Simulation
# This is the first simulation of two for estimating aerorol radiative effect

readonly COMPSET="F2010SC5-P3"   
readonly sstplus4K=false                     # true if to run plus4k experiment, otherwise false
readonly RESOLUTION="ne30pg2_r05_oECv3"      # or ne120pg2_r0125_oRRS18to6v3
readonly DESCRIPTOR="F2010-CTL.ne30pg2"      # This will be the main part of the casename

readonly CASE_GROUP="NGD.Convection"

# Code and compilation
readonly CHECKOUT="20210707"                       # Provide a timestamp for distinction
readonly BRANCH="wlin/atm/p3update_amipcompsets"   # update based off jacobshpundpnnl/atm/E3SMv2_alpha5_59_wP3v4
readonly DEBUG_COMPILE=false

# Run options
readonly MODEL_START_TYPE="initial"  # initial, continue
readonly START_DATE="2009-01-01"

# Case name
readonly CASE_NAME=${CHECKOUT}.${DESCRIPTOR}.${MACHINE}

# Set paths
readonly CODE_ROOT="/global/cscratch1/sd/wlin/share/E3SM.P3"               # where it contains 'components','cime', etc.
readonly CASE_ROOT="/global/cscratch1/sd/wlin/E3SM_testings/${CASE_NAME}"  

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'S_2x5_ndays', 'M_1x10_ndays', 'M80_1x10_ndays'
#  or 'production' for full simulation
#readonly run='M_1x2_nmonths'
readonly run='production'
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
  echo STOP_OPTION=${units}
  readonly STOP_N=${length}
  echo STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  echo RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false

else

  # Production simulation
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="M"
  readonly WALLTIME="24:00:00"
  readonly STOP_OPTION="nmonths"
  readonly STOP_N="25"
  readonly REST_OPTION="nmonths"
  readonly REST_N="3"
  readonly RESUBMIT="0"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history 
readonly HIST_OPTION="nmonths"
readonly HIST_N="1"

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
 cosp_lite = .true.

 avgflag_pertape = 'A','I'
 nhtfrq = 0,-6
 mfilt = 1,4
 fincl2 = 'PS','U','V','T','Q'

 inithist = 'MONTHLY'

 !--- aerosol forcing diagnostics !
 rad_diag_1 = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2', 'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4', 'N:CFC11:CFC11', 'N:CFC12:CFC12',

EOF

cat << EOF >> user_nl_elm

 finidat = '/global/cfs/cdirs/e3sm/shpundk/elm_inputdata/20201103.IELM.r05_oECv3.elm.r.0030-01-01-00000.nc'

EOF

}

patch_mpas_streams() {

echo
echo 'Modifying MPAS streams files' 
pushd ../run

# change steams.ocean file
patch streams.ocean << EOF
--- streams.ocean.00	2021-03-24 20:48:04.236324000 -0500
+++ streams.ocean	2021-03-24 20:54:34.145396766 -0500
@@ -533,7 +533,7 @@
 
 <stream name="timeSeriesStatsMonthlyOutput"
         type="output"
-        precision="single"
+        precision="double"
         io_type="pnetcdf"
         filename_template="20210305.v2beta3.piControl.ne30pg2_EC30to60E2r2.chrysalis.mpaso.hist.am.timeSeriesStatsMonthly.$Y-$M-$D.nc"
         filename_interval="00-01-00_00:00:00"
@@ -595,6 +595,7 @@
     <var name="longWaveHeatFluxDown"/>
     <var name="seaIceHeatFlux"/>
     <var name="shortWaveHeatFlux"/>
+    <var name="frazilTemperatureTendency"/>
     <var name="evaporationFlux"/>
     <var name="seaIceSalinityFlux"/>
     <var name="seaIceFreshWaterFlux"/>
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

    popd
}

#-----------------------------------------------------
create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
        echo $'\n----- Skipping create_newcase -----\n'
        return
    fi

    echo $'\n----- Starting create_newcase -----\n'

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

    # Reset CCSM_CO2_PPMV value to the one for 2010 CMIP6
    ./xmlchange CCSM_CO2_PPMV="388.717"

    # Conditional updates based on simulation config: sstplus4K or lustre stripe size,m pe-layout

    if [ "$sstplus4K" == "true" ]; then
     ./xmlchange SSTICE_DATA_FILENAME="$input_data_dir/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213_plus4K.nc"
    fi

    if [[ $RESOLUTION =~ "ne120" ]]; then
     # lfs setstripe -S 1m -c 64 `./xmlquery --value RUNDIR`
       ./xmlchange MAX_MPITASKS_PER_NODE="64"
       ./xmlchange NTASKS_WAV="32"
       ./xmlchange NTASKS_GLC="32"
       for comp in ATM LND ROF ICE OCN CPL; do
           ./xmlchange NTASKS_$comp="21600"
           ./xmlchange NTHRDS_$comp="2"
       done
    fi

    # Custom user_nl
    user_nl

    # reset pe-layout
    #cp ~/E3SM/Cases/prod/F-16nodes-chrys.env_mach_pes.xml env_mach_pes.xml

    # Finally, run CIME case.setup
    ./case.setup --reset

    # setstripe can only be done after case.setup (after RUNDIR is created)

    if [[ $MACHINE =~ "cori" ]] && [[ $RESOLUTION =~ "ne120" ]]; then
       lfs setstripe -S 1m -c 64 `./xmlquery --value RUNDIR`
    fi


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

    # TO DO: implement 'branch', 'hybrid'

    else
        echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
        exit 380
    fi

    # Patch mpas streams files. Which code base need it? Not needed anyway for F-case
    # patch_mpas_streams

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

