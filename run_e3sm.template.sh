#!/bin/bash -fe

# E3SM Coupled Model Group run_e3sm script template.
#
# Bash coding style inspired by:
# https://web.archive.org/web/20200620202413/http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

main() {

# For debugging, uncomment line below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=pm-gpu
# NOTE: The command below will return your default project on SLURM-based systems. 
# If you are not using SLURM or need a different project, remove the command and set it directly
readonly PROJECT="$(sacctmgr show user $USER format=DefaultAccount | tail -n1 | tr -d ' ')"

# Simulation
readonly COMPSET="F2010-SCREAMv1"
readonly RESOLUTION="ne30pg2_ne30pg2"
# BEFORE RUNNING : CHANGE the following CASE_NAME to desired value
readonly CASE_NAME="my_fancy_run.${COMPSET}.${RESOLUTION}.${MACHINE}"  
# If this is part of a simulation campaign, ask your group lead about using a case_group label
# readonly CASE_GROUP=""

# Code and compilation
# BEFORE RUNNING: CHANGE CHECKOUT to date string like 20240301
readonly CHECKOUT="latest"
readonly BRANCH="master"
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="0001-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
#readonly RUN_REFDIR=""
#readonly RUN_REFCASE=""
#readonly RUN_REFDATE=""   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

# Set paths
readonly CASE_ROOT="${PSCRATCH}/EAMxx/${CASE_NAME}"
readonly CODE_ROOT="/pscratch/sd/m/mahf708/e3sm-repo/test-pr"

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays',
#               'M_1x10_ndays', 'M2_1x10_ndays', 'M80_1x10_ndays', 'L_1x10_ndays'
#               * can replace XS, M, etc. with custom-XY with XY being the node count
#  or 'production' for full simulation
readonly run='4x1_1x6_ndays'
if [ "${run}" != "production" ]; then
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

# increase SCREAM_NUM_TRACERS from 10 to 11 (note that scream supports both 128 and 72 levels for low-res)
./xmlchange SCREAM_CMAKE_OPTIONS="SCREAM_NP 4 SCREAM_NUM_VERTICAL_LEV 128 SCREAM_NUM_TRACERS 11"

# add pompei to the list of aerosol processes
./atmchange mac_aero_mic::atm_procs_list+=pompei
# if you want to change the eruption date, you can do so with the following command
./atmchange atmosphere_processes::physics::mac_aero_mic::pompei::eruption_date="0001-01-02-00000"
# if you want to change the eruption radius, you can do so with the following command
./atmchange atmosphere_processes::physics::mac_aero_mic::pompei::plume_radius_in_km=1000.0

# create yaml file (or save it elsewhere)
cat << EOF >> output_file.yaml
%YAML 1.1
---
filename_prefix: tutorial_output.eamxx.h
Averaging Type: Average
Max Snapshots Per File: 1
track_fill: true
Fields:
  Physics PG2:
    Field Names:
    # 3D vars
    - ash
    - T_mid
    - qv
    - RelativeHumidity
    - qc
    - qi
    - qr
    - qm
    - nc
    - ni
    - nr
    - cldfrac_tot_for_analysis
    - cldfrac_ice_for_analysis
    - cldfrac_liq
    - omega
    - U
    - V
    - z_mid
    - p_mid
    - tke
    # 2D vars
    - ash_column
    - ash_at_model_top
    - ash_at_model_bot
    - ash_at_400hPa
    - ash_at_700hPa
    - ash_at_800hPa
    - ash_at_900hPa
    - ash_at_990hPa
    - SW_flux_up_at_model_top
    - SW_flux_dn_at_model_top
    - LW_flux_up_at_model_top
    - SW_clrsky_flux_up_at_model_top
    - LW_clrsky_flux_up_at_model_top
    - SW_flux_dn_at_model_bot
    - SW_clrsky_flux_dn_at_model_bot
    - SW_flux_up_at_model_bot
    - SW_clrsky_flux_up_at_model_bot
    - LW_flux_dn_at_model_bot
    - LW_clrsky_flux_dn_at_model_bot
    - LW_flux_up_at_model_bot
    - LongwaveCloudForcing
    - ShortwaveCloudForcing
    - ps
    - SeaLevelPressure
    - T_2m
    - qv_2m
    - surf_radiative_T
    - VapWaterPath
    - IceWaterPath
    - LiqWaterPath
    - RainWaterPath
    - ZonalVapFlux
    - MeridionalVapFlux
    - surf_evap
    - surf_sens_flux
    - surface_upward_latent_heat_flux
    - precip_liq_surf_mass_flux
    - precip_ice_surf_mass_flux
    - landfrac
    - ocnfrac
    - PotentialTemperature_at_700hPa
    - PotentialTemperature_at_1000hPa
    - omega_at_500hPa
    - RelativeHumidity_at_700hPa
output_control:
  Frequency: 1
  frequency_units: ndays
  MPI Ranks in Filename: false
EOF

# add the output file to the list of output streams
./atmchange output_yaml_files="./output_file.yaml"

cat << EOF >> user_nl_elm
finidat = ''
hist_empty_htapes=.true.
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
			--pecount ${layout}
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
			--pecount ${layout}
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

    # Build with COSP, except for a data atmosphere (datm) or "scream"
    if [ `./xmlquery --value COMP_ATM` == "datm" ] || [ `./xmlquery --value COMP_ATM` == "scream" ]; then
      echo $'\nThe specified configuration uses a data atmosphere or SCREAM, so will not add COSP to CAM_CONFIG_OPTS'
      echo $'If you want to use the COSP simulator in EAMxx, you can add it via atmchange as a process.\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # Finally, run CIME case.setup
    ./case.setup --reset

    # Custom user_nl (if we use atmchange inside user_nl, we need to call it after case.setup?)
    user_nl

    popd
}

#-----------------------------------------------------
custom_pelayout() {

if [[ ${PELAYOUT} == custom-* ]];
then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    # Number of cores per node (machine specific)
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
	 echo $'\n----- Preview namelists -----\n'
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
    local this_script_name=$( basename -- "$0"; )
    local this_script_dir=$( dirname -- "$0"; )
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp "${this_script_dir}/${this_script_name}" ${script_provenance_dir}/${script_provenance_name}

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

