#!/bin/bash -fe
# E3SM+GCAM v3 SSP245 production run script for chrysalis
# The only thing that must be set correctly below by the user is the MACHINE variable
# If you alread have a code clone, set CLONE_NAME and CODE_PARENT appropriately below
# If you want to fetch the code then set and uncomment 'code and compilation' below,
#    and set do_fetch_code=true
# You may consider setting a more specific CASE_NAME below
# If desired, the user can change other things below, but this script will perform as desired

main() {

# For debugging, uncomment libe below
#set -x

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_modify_pe_layout=true
do_case_setup=true
do_case_build=true
do_case_submit=true

# --- Configuration flags ----

# Machine and project
# Be sure to set your machine here!
readonly PROJECT="e3sm"
#readonly MACHINE="chrysalis"
#readonly MACHINE="compy"
readonly MACHINE="pm-cpu"

# if true this increases the active component ntasks to ~5400 from ~3200
ehc_pe_xl=false

# set the machine inputdata directory
if [ "${MACHINE}" == "chrysalis" ]; then
   readonly din_loc_root=/lcrc/group/e3sm/data/inputdata
fi
if [ "${MACHINE}" == "compy" ]; then
   readonly din_loc_root=/compyfs/inputdata
fi
if [ "${MACHINE}" == "pm-cpu" ]; then
   din_loc_root=/global/cfs/cdirs/e3sm/inputdata
fi

# Simulation
readonly MYDATE=$(date '+%Y%m%d%H') # use current date if MYDATE is not set to a specific date
# export COMPSET=SSP245_EAM%CMIP6_ELM%TOPCNPRDCTCBCPHS_MPASSI%PRES_DOCN%DOM_SROF_SGLC_SWAV_GCAM_BGC%LNDATM
readonly COMPSET="SSP245_ZATM_BGC" # see long name above
readonly RESOLUTION="ne30pg2_f09_oEC60to30v3" 
readonly CASE_NAME="${COMPSET}_${RESOLUTION}_${MYDATE}_v3_pr_test"
# readonly CASE_GROUP="E3SM_GCAM"

# Code and compilation
#readonly BRANCH="master" 
#readonly giac_branch="master"
#readonly gcam_branch="e3sm-integration"

# compile with debug?
readonly DEBUG_COMPILE=FALSE

# Run options
readonly MODEL_START_TYPE="hybrid"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="2015-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
readonly RUN_REFCASE="20260303_I20TREAMELMCNPRDCTCBCPHSBGC_${RESOLUTION}" 
readonly RUN_REFDATE="2015-01-01"
readonly RUN_REFDIR="$din_loc_root/e3sm_init/${RUN_REFCASE}/${RUN_REFDATE}-00000"
readonly MPASSI_CONFIG_START=2015-01-01_0

# Set paths
# note that fetch_code below will use CLONE_NAME in place of E3SM, so that CODE_ROOT points to the root of the cloned repository, and not to a parent directory containing multiple clones
readonly CLONE_DATE=$(date '+%d%B%Y')

#readonly CLONE_NAME="e3sm_gcam_${CLONE_DATE}"
export CLONE_NAME="e3sm_gcam_april2026_for_rebase"

readonly CODE_PARENT="${HOME}/e3sm"
readonly CODE_ROOT="${CODE_PARENT}/${CLONE_NAME}"

# set the machine scratch directory
if [ "${MACHINE}" == "chrysalis" ]; then
   readonly CASE_ROOT="/lcrc/group/e3sm/$USER/e3sm_scratch/${CASE_NAME}"
fi
if [ "${MACHINE}" == "compy" ]; then
   readonly CASE_ROOT="/compyfs/${USER}/e3sm_scratch/${CASE_NAME}"
fi
if [ "${MACHINE}" == "pm-cpu" ]; then
   readonly CASE_ROOT="${SCRATCH}/e3sm_scratch/${CASE_NAME}"
fi

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays', 
#               'M_1x10_ndays', 'M2_1x10_ndays', 'M80_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation
# readonly run='XS_2x5_ndays'
# readonly run='XL_2x12_nmonths' 
# readonly run='L_1x12_nmonths'
# readonly run='L_1x5_nyears'
# readonly run='M_12x5_nmonths'
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
  readonly PELAYOUT="EHC-ELM-EAM"
  readonly WALLTIME="14:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="6"
  readonly REST_OPTION="nyears"
  readonly REST_N="1"
  readonly RESUBMIT="0"
  readonly DO_SHORT_TERM_ARCHIVING=false
  readonly IFDEBUG="regular"
fi

# Coupler history 
readonly HIST_OPTION="nyears"
readonly HIST_N="1"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

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

ncd_string="${RUN_REFDIR}/${RUN_REFCASE}.eam.i.${RUN_REFDATE}-00000.nc"

cat << EOF >> user_nl_eam
 co2_print_diags_timestep               = .true.
 co2_print_diags_monthly                = .true.
 co2_print_diags_total                  = .true.
 cflx_cpl_opt=1 
 
 ncdata		= '${ncd_string}'
EOF

# don't need finidat because doing hybrid run from initial conditions refcase restart
cat << EOF >> user_nl_elm

 hist_mfilt = 1, 365, 1
 hist_nhtfrq = 0, -24, 0
 hist_dov2xy = .true., .true., .false.
 hist_fincl2 = 'TBOT', 'TREFMXAV', 'TREFMNAV', 'RAIN', 'SNOW', 'SNOWDP'
 hist_fincl3 = 'GPP', 'ER', 'HR', 'NPP'
EOF

cat << EOF >> user_nl_gcam
!read_scalars = .true.
!scalar_source_dir = ''
EOF

cat << EOF >> user_nl_cpl
EOF

cat >> user_nl_mpassi << EOF
&seaice_model
 config_start_time = '$MPASSI_CONFIG_START'
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
    # recall that CODE_ROOT is set to the root of the cloned repository, not the parent directory of the clone
    local path=${CODE_ROOT}
    local repo=e3sm

    echo "Cloning $repo repository branch $BRANCH into $path"
    if [ -d "${path}" ]; then
        echo "ERROR: Directory already exists. Not overwriting"
        exit 20
    fi

    # This will put repository E3SM into CLONE_NAME, within the CODE_PARENT directory
    # For example, if CODE_PARENT=/home/user/e3sm and CLONE_NAME=clone1, then this will put the code in /home/user/e3sm/clone1
    git clone git@github.com:E3SM-Project/${repo}.git ${CLONE_NAME}
    cd ${CODE_ROOT}
    
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
    # to reduce clone size, exclude submodule history by adding --depth=1
    git submodule update --init --recursive

    # Check out submodule branches if necessary
    #cd components/gcam/src
    #git checkout $giac_branch
    #cd iac/gcam
    #git checkout $gcam_branch
    #cd ../../../../..

    cd ..

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

    # machine specific settings
    if [ `./xmlquery --value MACH` == chrysalis ]; then
        ./xmlchange MAX_TASKS_PER_NODE=64
        ./xmlchange MAX_MPITASKS_PER_NODE=64
        if [ "${ehc_pe_xl}" == "true" ]; then
           readonly nnodes=85
        else
           readonly nnodes=52
        fi
    fi
    if [ `./xmlquery --value MACH` == compy ]; then
        if [ "${ehc_pe_xl}" == "true" ]; then
           readonly nnodes=136
        else
           readonly nnodes=80
        fi
	    ./xmlchange MAX_TASKS_PER_NODE=40
        ./xmlchange MAX_MPITASKS_PER_NODE=40
    fi
    if [ `./xmlquery --value MACH` == pm-cpu ]; then
        if [ "${ehc_pe_xl}" == "true" ]; then
           readonly nnodes=43
        else
           readonly nnodes=26
        fi
        ./xmlchange MAX_TASKS_PER_NODE=128
        ./xmlchange MAX_MPITASKS_PER_NODE=128
    fi

    readonly ppn=`./xmlquery MAX_MPITASKS_PER_NODE --value`

    ./xmlchange NTASKS_WAV=1
    ./xmlchange NTASKS_GLC=1
    ./xmlchange NTASKS_ESP=1
    ./xmlchange NTASKS_IAC=1

    ./xmlchange ROOTPE=0
    ./xmlchange NTHRDS=1

    ./xmlchange NTHRDS_ATM=1
    ./xmlchange NTHRDS_CPL=1
    ./xmlchange NTHRDS_OCN=1
    ./xmlchange NTHRDS_WAV=1
    ./xmlchange NTHRDS_GLC=1
    ./xmlchange NTHRDS_ICE=1
    ./xmlchange NTHRDS_ROF=1
    ./xmlchange NTHRDS_LND=1
    ./xmlchange NTHRDS_IAC=1

    ./xmlchange ROOTPE_ATM=0
    ./xmlchange ROOTPE_CPL=0
    ./xmlchange ROOTPE_OCN=0
    ./xmlchange ROOTPE_WAV=0
    ./xmlchange ROOTPE_GLC=0
    ./xmlchange ROOTPE_ICE=0
    ./xmlchange ROOTPE_ROF=0
    ./xmlchange ROOTPE_LND=0
    ./xmlchange ROOTPE_IAC=0

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
        #  ./xmlchange NTASKS_OCN=5440
        #  ./xmlchange NTASKS_ICE=5440
         ./xmlchange NTASKS_OCN=5400
         ./xmlchange NTASKS_ICE=5400
         ./xmlchange NTASKS_ROF=5440
         ./xmlchange NTASKS_LND=5440
    elif [[ ${PELAYOUT} == "EHC-ELM-EAM" ]]; then
         echo $'\n----- changing PE layout to EHC-ELM-EAM -----\n'
         ./xmlchange NTASKS_ATM=$(($ppn * $nnodes))
         ./xmlchange NTASKS_CPL=$(($ppn * $nnodes))
         ./xmlchange NTASKS_LND=$(($ppn * $nnodes))
         ./xmlchange NTASKS_ICE=256
         ./xmlchange NTASKS_OCN=256
         ./xmlchange NTASKS_ROF=1
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

    # Build with COSP, except for a data atmosphere (datm)
    if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then 
      echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    readonly input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # Custom user_nl
    user_nl

    # these are required here for the case to be set up properly, even if some are runtime settings

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

    else
        echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
        exit 380
    fi

    ./xmlchange SAVE_TIMING=TRUE

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

        # clean the build if an executable already exists, just to be safe
        if [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
            ./case.build --clean-all
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

 

    # Patch mpas streams files
    patch_mpas_streams

    # # modify to use priority queue for simulation campaigns
    # if [ "${run}" = "production" ]; then
    #     ./xmlchange --file env_workflow.xml --id CHARGE_ACCOUNT --val "priority"
    #     ./xmlchange --file env_workflow.xml --id PROJECT --val "priority"
    #     ./xmlchange --file env_workflow.xml --id JOB_QUEUE --force --val "priority"
    # fi

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
    local this_script_name=$(realpath "$0")
    local script_provenance_name=$(basename ${this_script_name}).`date +%Y%m%d-%H%M%S`
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
