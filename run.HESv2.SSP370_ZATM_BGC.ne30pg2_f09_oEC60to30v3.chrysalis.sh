#!/bin/bash -fe
# E3SM+GCAM SSP370 production run script for chrysalis

main() {

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=chrysalis
readonly PROJECT="e3sm"

# Simulation
readonly MYDATE=$(date '+%Y%m%d%H') # use current date if MYDATE is not set
#export COMPSET=SSP370_EAM%CMIP6_ELM%CNPRDCTCBC_MPASSI%PRES_DOCN%DOM_SROF_SGLC_SWAV_GCAM_BGC%LNDATM
export COMPSET="SSP370_ZATM_BGC"
readonly RESOLUTION="ne30pg2_f09_oEC60to30v3" # non-default grids are: atm:ne30np4.pg2  lnd:0.9x1.25  ocnice:oEC60to30v3  rof:null  glc:null  wav:null   mask is: oEC60to30v3
readonly CASE_NAME="${COMPSET}_${RESOLUTION}_${MYDATE}"
# readonly CASE_GROUP="E3SM_GCAM"

# Code and compilation
#readonly CHECKOUT="${MYDATE}"
#readonly BRANCH="main" 
#$(git rev-parse HEAD) #"9f0094d40f039769ba68ed4e89810516838ad43f"
#readonly CHERRY=( "7e4d1c9fec40ce1cf2c272d671f5d9111fa4dea7" "a5b1d42d7cd24924d0dbda95e24ad8d4556d93f1" ) # PR4349
readonly DEBUG_COMPILE=TRUE

# Run options
readonly MODEL_START_TYPE="hybrid"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="2015-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
readonly RUN_REFCASE="20241204_I20TREAMELMCNPRDCTCBCBGC_${RESOLUTION}" 
readonly RUN_REFDATE="2015-01-01"
readonly RUN_REFDIR="/lcrc/group/e3sm/ac.eva.sinha/E3SM_GCAM_lnd_init/${RUN_REFCASE}"

# Set paths
readonly CODE_ROOT="${HOME}/code/E3SM_GCAM/E3SM"
readonly CASE_ROOT="/lcrc/globalscratch/ac.sfeng1/${CASE_NAME}"

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays', 
#               'M_1x10_ndays', 'M2_1x10_ndays', 'M80_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation
readonly run='XS_2x5_ndays'
#readonly run='XL_12x2_nmonths' 
#readonly run='L_12x2_nmonths'
#readonly run='M_12x5_nmonths'

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
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false
  
  readonly WALLTIME="4:00:00"
  readonly IFDEBUG="debug"
else

  # Production simulation
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="XL"
  readonly WALLTIME="30:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="40"
  readonly REST_OPTION="nyears"
  readonly REST_N="5"
  readonly RESUBMIT="2"
  readonly DO_SHORT_TERM_ARCHIVING=false
  readonly IFDEBUG="compute"
fi

# Coupler history 
readonly HIST_OPTION="nyears"
readonly HIST_N="5"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_modify_pe_layout=true
do_case_setup=true
do_case_build=true
# do_case_submit=true
do_case_submit=false
# --- Now, do the work ---

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

#change PE layout
modify_pe_layout

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
 co2_conserv_error_tol_per_year = 1.e-5
 
 ncdata = '20241204_I20TREAMELMCNPRDCTCBCBGC_ne30pg2_f09_oEC60to30v3.eam.i.2015-01-01-00000.nc'
 
 hist_mfilt = 1
 hist_nhtfrq = -24
EOF

# Setting do_harvest == .false. because the iac takes care of this
cat << EOF >> user_nl_elm
 do_budgets = .true.
 do_harvest = .false.

 finidat = '20241204_I20TREAMELMCNPRDCTCBCBGC_ne30pg2_f09_oEC60to30v3.elm.r.2015-01-01-00000.nc'
 
 hist_dov2xy = .true.
 
 hist_mfilt = 1
 hist_nhtfrq = -24
EOF

cat << EOF >> user_nl_gcam

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
	    --queue ${IFDEBUG} \
        --pecount ${PELAYOUT}

    if [ $? != 0 ]; then
      echo $'\nNote: if create_newcase failed because sub-directory already exists:'
      echo $'  * delete old case_script sub-directory'
      echo $'  * or set do_newcase=false\n'
      exit 35
    fi
}

modify_pe_layout() {

    if [ "${do_modify_pe_layout,,}" != "true" ]; then
        echo $'\n----- Skipping changing PE-layout -----\n'
        return
    fi

    pushd ${CASE_SCRIPTS_DIR}

    ./xmlchange MAX_MPITASKS_PER_NODE=64
    ./xmlchange MAX_TASKS_PER_NODE=64

    ./xmlchange NTASKS_WAV=1
    ./xmlchange NTASKS_GLC=1

    ./xmlchange NTHRDS_ATM=1
    ./xmlchange NTHRDS_CPL=1
    ./xmlchange NTHRDS_OCN=1
    ./xmlchange NTHRDS_WAV=1
    ./xmlchange NTHRDS_GLC=1
    ./xmlchange NTHRDS_ICE=1
    ./xmlchange NTHRDS_ROF=1
    ./xmlchange NTHRDS_LND=1

    ./xmlchange ROOTPE_ATM=0
    ./xmlchange ROOTPE_CPL=0
    ./xmlchange ROOTPE_OCN=0
    ./xmlchange ROOTPE_WAV=0
    ./xmlchange ROOTPE_GLC=0
    ./xmlchange ROOTPE_ICE=0
    ./xmlchange ROOTPE_ROF=0
    ./xmlchange ROOTPE_LND=0

    if [[ ${PELAYOUT} == "XS" ]]; then
        echo $'\n----- changing PE layout to XS -----\n'
	     ./xmlchange NTASKS_ATM=320
        ./xmlchange NTASKS_CPL=320
        ./xmlchange NTASKS_OCN=320
        ./xmlchange NTASKS_ICE=256
        ./xmlchange NTASKS_WAV=256
        ./xmlchange NTASKS_ROF=256
        ./xmlchange NTASKS_LND=320    
    elif [[ ${PELAYOUT} == "S" ]]; then
        echo $'\n----- changing PE layout to S -----\n'
        ./xmlchange NTASKS_ATM=1080
	      ./xmlchange NTASKS_CPL=1280
        ./xmlchange NTASKS_OCN=1280
        ./xmlchange NTASKS_ICE=1280
        ./xmlchange NTASKS_ROF=1280
        ./xmlchange NTASKS_LND=1280
    elif [[ ${PELAYOUT} == "M" ]]; then
        echo $'\n----- changing PE layout to M -----\n'
	      ./xmlchange NTASKS_ATM=1350
        ./xmlchange NTASKS_CPL=1408
        ./xmlchange NTASKS_OCN=1408
        ./xmlchange NTASKS_ICE=1408
        ./xmlchange NTASKS_ROF=1408
        ./xmlchange NTASKS_LND=1408
    elif [[ ${PELAYOUT} == "L" ]]; then
        echo $'\n----- changing PE layout to L -----\n'
        ./xmlchange NTASKS_ATM=2700
	      ./xmlchange NTASKS_CPL=2752
        ./xmlchange NTASKS_OCN=2752
        ./xmlchange NTASKS_ICE=2752
        ./xmlchange NTASKS_ROF=2752
        ./xmlchange NTASKS_LND=2752
     elif [[ ${PELAYOUT} == "XL" ]]; then
         echo $'\n----- changing PE layout to XL -----\n'
         ./xmlchange NTASKS_ATM=5400
         ./xmlchange NTASKS_CPL=5440
         ./xmlchange NTASKS_OCN=5440
         ./xmlchange NTASKS_ICE=5440
         ./xmlchange NTASKS_ROF=5440
         ./xmlchange NTASKS_LND=5440
    else
        echo 'ERROR: $PELAYOUT = '${PELAYOUT}' but no layout exists for this setting.'
        exit 297
    fi

    popd
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
    
    # using custom layout for cori-kn
    ./xmlchange --file env_mach_pes.xml --id MAX_TASKS_PER_NODE --val 64
    ./xmlchange --file env_mach_pes.xml --id MAX_MPITASKS_PER_NODE --val 64

    declare -a comps=("ATM" "CPL" "LND" "ICE" "OCN" "GLC" "ROF" "WAV" "ESP" )
    for thiscomp in "${comps[@]}"; do
      ./xmlchange --file env_mach_pes.xml --id NTHRDS_$thiscomp --val 1
    done

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

    # ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR} # may not be neccessary 
    # docn setup
    ./xmlchange --id PIO_TYPENAME  --val "pnetcdf"
    ./xmlchange SSTICE_DATA_FILENAME=${input_data_dir}/ocn/docn7/SSTDATA/sst_ice_GFDL-ESM4_ssp245_r2i1p1f1_gr_201501-210012_land_interpolated.nc
    ./xmlchange SSTICE_YEAR_START=2015
    ./xmlchange SSTICE_YEAR_END=2100
    ./xmlchange SSTICE_YEAR_ALIGN=2015

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

    else

        echo $'\n----- Starting case_build -----\n'
        echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"

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
        ./xmlchange RUN_TYPE=${MODEL_START_TYPE}
        ./xmlchange GET_REFCASE=${GET_REFCASE}
	    ./xmlchange RUN_REFDIR=${RUN_REFDIR}
        ./xmlchange RUN_REFCASE=${RUN_REFCASE}
        ./xmlchange RUN_REFDATE=${RUN_REFDATE}
        echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE} 
        echo '$RUN_REFDIR = '${RUN_REFDIR}
        echo '$RUN_REFCASE = '${RUN_REFCASE}
        echo '$RUN_REFDATE = '${START_DATE}

        ln -sf  ${RUN_REFDIR}/${RUN_REFCASE}.eam.i.${RUN_REFDATE}-00000.nc ${CASE_RUN_DIR}/.
        ln -sf  ${RUN_REFDIR}/${RUN_REFCASE}.elm.r.${RUN_REFDATE}-00000.nc ${CASE_RUN_DIR}/.

        else
        echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
        exit 380
    fi

    # Patch mpas streams files
    patch_mpas_streams

    # modify to use priority queue 
    if [ "${run}" = "production" ]; then
        ./xmlchange --file env_workflow.xml --id CHARGE_ACCOUNT --val "priority"
        ./xmlchange --file env_workflow.xml --id PROJECT --val "priority"
        ./xmlchange --file env_workflow.xml --id JOB_QUEUE --force --val "priority"
    fi

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