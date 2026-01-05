#!/bin/bash -fe

# E3SM Coupled Model Group run_e3sm script template.
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
readonly MACHINE=pm-cpu
# BEFORE RUNNING:  CHANGE this to your project
readonly PROJECT="m4942"
readonly QUEUE="debug"  # 'regular', 'debug', 'premium'
readonly COMPILER="gnu"

# Simulation
readonly COMPSET="F2010"
# Full ne30: "ne30pg2_r05_IcoswISC30E3r5" Test ne4: "ne4pg2_oQU480"
readonly RESOLUTION="ne4pg2_oQU480" 
#readonly RESOLUTION="ne30pg2_r05_IcoswISC30E3r5"
# BEFORE RUNNING : CHANGE the following CASE_NAME to desired value
readonly CASE_NAME="e3sm2025_mlmicro_dnn_emulator_withinputfilters5"
# If this is part of a simulation campaign, ask your group lead about using a case_group label
# readonly CASE_GROUP=""

# Code and compilation
# BEFORE RUNNING: CHANGE CHECKOUT to date string like 20240301
readonly CHECKOUT="20240506"
readonly BRANCH="master"
readonly CHERRY=( )
readonly DEBUG_COMPILE=FALSE

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="2001-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=FALSE
#readonly RUN_REFDIR=""
#readonly RUN_REFCASE=""
#readonly RUN_REFDATE=""   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

# Set paths
# Use local clone (read/write) to avoid permission issues
readonly CODE_ROOT="${HOME}/e3sm_ftorch"
#readonly CODE_ROOT="/global/cfs/cdirs/m4549/code/e3sm2025_ftorch_mlmicro"  
#readonly CODE_ROOT="/pscratch/sd/p/plma/shared/for_andrew/E3SM"
#readonly CASE_ROOT="/pscratch/sd/a/agett/e3sm_scratch/${CASE_NAME}"
readonly CASE_ROOT="/pscratch/sd/d/dvpatel/mlmicrophysics_project/e3sm_ftorch_scratch/${CASE_NAME}"
rm -rf $CASE_ROOT

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays',
#               'M_1x10_ndays', 'M2_1x10_ndays', 'M80_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation
readonly run="XS_1x1_ndays"
#readonly run="production"

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
  readonly WALLTIME="0:30:00"
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
  readonly PELAYOUT="custom-43"
  readonly WALLTIME="06:00:00"
  readonly STOP_OPTION="nmonths"
  readonly STOP_N="12"
  readonly REST_OPTION="nmonths"
  readonly REST_N="12"
  readonly RESUBMIT="1"
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

#           'P3_mu_c','P3_lamc','P3_nr','P3_lamr',

# To add a new emulator: replace p3_torchscript_warm_rain_emulator_file with the new file.
# To activate it (right now set for training running tau, set p3_warm_rain_method = 'ftorch_emulator')

user_nl() {

cat << EOF >> user_nl_eam

 p3_warm_rain_method = 'ftorch_emulator'

 p3_torchscript_warm_rain_emulator_file = '/pscratch/sd/d/dvpatel/mlmicrophysics_project/trained_emulator_files/emulator107666_torchscript_qctauin1e-6_cloud1e-2.pt'
 p3_tau_kernel_file = '/global/cfs/cdirs/m4942/e3sm/emulators/old_emulator_files/KBARF_tau_kernel.dat' 
 p3_stochastic_emulated_filename_quantile = '/global/cfs/cdirs/m4942/e3sm/emulators/old_emulator_files/quantile_neural_net_fortran.nc'
 p3_stochastic_emulated_filename_input_scale = '/global/cfs/cdirs/m4942/e3sm/emulators/old_emulator_files/input_quantile_scaler.nc'
 p3_stochastic_emulated_filename_output_scale = '/global/cfs/cdirs/m4942/e3sm/emulators/old_emulator_files/output_quantile_scaler.nc'

 cosp_lite = .false.

 empty_htapes = .true.

 ! -- chemUCI settings ------------------
 history_chemdyg_summary = .true.
 history_gaschmbudget_2D = .false.
 history_gaschmbudget_2D_levels = .false.
 history_gaschmbudget_num = 6 !! no impact if  history_gaschmbudget_2D = .false.

 ! -- MAM5 settings ------------------
 is_output_interactive_volc = .true.

 nhtfrq = 0,-6
 mfilt  = 1,12
 avgflag_pertape = 'A','I'
 fincl1 = 'T','RELHUM','FLNT','FSNT','SWCF','LWCF','PRECT','PRECL','CLDTOT',
          'TGCLDLWP','TGCLDIWP','RAINQM','CLDLIQ','NUMLIQ','NUMRAI','CLDICE','CLOUD',
          'P3_qc2qr_accret_tend','P3_qc2qr_autoconv_tend', 'P3_nc_accret_tend', 
          'P3_nc2nr_autoconv_tend','P3_nc_selfcollect_tend', 'P3_nr_selfcollect_tend', 
          'P3_qctend_TAU','P3_nctend_TAU','P3_qrtend_TAU','P3_nrtend_TAU',
          'P3_qctend_TAU_raw','P3_nctend_TAU_raw','P3_qrtend_TAU_raw','P3_nrtend_TAU_raw'
 fincl2 = 'RELHUM','PRECT','PRECL','CLDTOT',
          'TGCLDLWP','TGCLDIWP','RAINQM','CLDLIQ','NUMLIQ','NUMRAI','CLDICE',
          'CLOUD','FREQR','RHO_CLUBB','T',
          'P3_mu_c','P3_lamc','P3_nr','P3_lamr',
          'P3_qctend_TAU_raw','P3_nctend_TAU_raw','P3_qrtend_TAU_raw','P3_nrtend_TAU_raw',
          'P3_qc_in_TAU','P3_nc_in_TAU','P3_qr_in_TAU','P3_nr_in_TAU',
          'P3_qc_out_TAU','P3_nc_out_TAU','P3_qr_out_TAU','P3_nr_out_TAU'

EOF
}

#finidat = '/pscratch/sd/p/plma/v3.LR.amip_0101_archive/rest/2001-01-01-00000/v3.LR.amip_0101.elm.r.2001-01-01-00000.nc'

#cat << EOF >> user_nl_elm
#finidat = '/global/cfs/cdirs/e3sm/inputdata/e3sm_init/20241011.v3.LR.HES-CBGC01-newIC.spinupIcepack/0051-01-01-00000/20241011.v3.LR.HES-CBGC01-newIC.spinupIcepack.elm.r.0051-01-01-00000.nc'
#hist_dov2xy = .true.,.true.
# check_finidat_year_consistency = .false.
# check_finidat_pct_consistency = .false.
# check_dynpft_consistency = .false.
# create_crop_landunit = .false.
#EOF

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
            --compiler ${COMPILER} \
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
            --compiler ${COMPILER} \
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
	 echo $'\n----- Preview namelists -----\n'
        ./preview_namelists

    fi

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Set Queue
    ./xmlchange JOB_QUEUE=${QUEUE}

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

# =====================================================
# Custom PE layout: custom-N where N is number of nodes
# =====================================================

custom_pelayout(){

if [[ ${PELAYOUT} == custom-* ]];
then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    # Number of cores per node (machine specific)
#    if [ "${MACHINE}" == "chrysalis" ]; then
        ncore=64
        hthrd=1  # hyper-threading
#    else
#        echo 'ERROR: MACHINE = '${MACHINE}' is not supported for current custom PE layout setting.'
#        exit 400
#    fi

    # Extract number of nodes
    tmp=($(echo ${PELAYOUT} | tr "-" " "))
    nnodes=${tmp[1]}

    # Applicable to all custom layouts
    pushd ${CASE_SCRIPTS_DIR}

    ./xmlchange NTASKS=$(( $nnodes * $ncore ))
    ./xmlchange NTHRDS=1
    ./xmlchange ROOTPE=0
    ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
    ./xmlchange MAX_TASKS_PER_NODE=$(( $ncore * $hthrd))

    popd

fi

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
