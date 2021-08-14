#!/bin/csh
date

set echo verbose


set fetch_code    = 0   # 0 = No, >0 = Yes
set compile_model = 1   # 0 = No, >0 = Yes
set run_model     = 1   # 0 = No, >0 = Yes

####################################################################
# Fetch code
####################################################################
setenv CODE_ROOT /compyfs/${user}/code
setenv CCSM_TAG  e3smv2_1d66ec33
setenv CCSM_ROOT ${CODE_ROOT}/${CCSM_TAG}

setenv BRANCH $CCSM_TAG 

if ($fetch_code > 0) then

   cd ${CODE_ROOT}
   git git@github.com:E3SM-Project/E3SM.git $CCSM_TAG

   # Setup git hooks
   rm -rf .git/hooks
   git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
   git config commit.template .git/hooks/commit.template
   git checkout ${BRANCH}
   # Bring in all submodule components
   git submodule update --init --recursive

endif

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv COMPSET    "F2010"
setenv RESOLUTION "ne30pg2_r05_oECv3"
setenv MACH       compy
setenv PROJECT    "esmd"

setenv CASE_NAME        RHI125_${MACH}_${COMPSET}_${RESOLUTION}_${CCSM_TAG}
setenv PTMP             /compyfs/$user/cond_diag_scratch

setenv CASE_ROOT        $PTMP/$CASE_NAME/cases
setenv CASE_RUN_DIR     $PTMP/$CASE_NAME/run
setenv CASE_BUILD_DIR   $PTMP/$CASE_NAME/build

####################################################################
# Run options setup 
####################################################################
set MODEL_START_TYPE = "initial"
set MODEL_START_DATE = "0001-01-01"

# Additional options for 'branch' and 'hybrid'
set GET_REFCASE      = FALSE
set RUN_REFDIR       = ""  # reference case run directory  
set RUN_REFCASE      = ""  # reference case name
set RUN_REFDATE      = ""  # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

set DEBUG            = 'FALSE'
set JOB_QUEUE        = "short"
set WALLTIME         = "2:00:00"
set STOP_OPTION      = "nmonths"
set STOP_N           = "1"
set RESUBMIT         = "0"
set REST_OPTION      = "nmonths"
set REST_N           = "1" 

set PELAYOUT         = "L"
set NTASKS           = 900
set NTHRDS           = 1 

#
####################################################################
# Compile model
####################################################################
if ($compile_model > 0 || $run_model > 0) then

    rm -rvf $CASE_ROOT

    ${CCSM_ROOT}/cime/scripts/create_newcase \
                    --case $CASE_NAME \
                    --output-root ${CASE_ROOT} \
                    --script-root ${CASE_ROOT} \
                    --handle-preexisting-dirs u \
                    --compset $COMPSET \
                    --res $RESOLUTION \
                    --machine $MACH \
                    --project ${PROJECT} \
                    --walltime ${WALLTIME} \
                    --pecount ${PELAYOUT}

#====================================================================
# set up case
#====================================================================

   ###./create_newcase -list grids

   cd $CASE_ROOT

   ./xmlchange -file env_run.xml      -id RUNDIR     -val $CASE_RUN_DIR
   ./xmlchange -file env_build.xml    -id EXEROOT    -val $CASE_BUILD_DIR

   # Model input directory 
   #./xmlchange -file env_run.xml      -id DIN_LOC_ROOT          -val $CSMDATA
   #./xmlchange -file env_run.xml      -id DIN_LOC_ROOT_CLMFORC  -val $CSMDATA

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $NTASKS
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $NTHRDS
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $NTASKS
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $NTHRDS
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $NTASKS
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $NTHRDS
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $NTASKS
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $NTHRDS
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $NTASKS
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $NTHRDS
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $NTASKS
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $NTHRDS
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $NTASKS
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $NTHRDS
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $NTASKS
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $NTHRDS
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val '0'

   ./case.setup --reset

#====================================================================
# Compile 
#====================================================================
   cd $CASE_ROOT

#  ./xmlchange -file env_build.xml -id DEBUG       -val 'TRUE'
#  ./xmlchange -file env_build.xml -id PIO_VERSION -val '1'


   if ($compile_model > 0) then # Build the model

      #rm -rvf $CASE_BUILD_DIR
      ./case.build

   else # mark as already built

      ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

   endif

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then
#------------------
## set environment
#------------------

cd $CASE_ROOT

./xmlchange  -file env_run.xml  -id  RUN_STARTDATE      -val $MODEL_START_DATE
./xmlchange  -file env_run.xml  -id  STOP_N             -val $STOP_N
./xmlchange  -file env_run.xml  -id  STOP_OPTION        -val $STOP_OPTION
./xmlchange  -file env_run.xml  -id  REST_N             -val $REST_N
./xmlchange  -file env_run.xml  -id  REST_OPTION        -val $REST_OPTION
./xmlchange  -file env_run.xml  -id  DOUT_S             -val 'FALSE'

./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val $WALLTIME
./xmlchange  -file env_batch.xml -id JOB_QUEUE          -val $JOB_QUEUE

if ( $RESUBMIT > 0 ) then
  ./xmlchange  -file env_run.xml  -id  RESUBMIT           -val $RESUBMIT
endif 
 
if ( $MODEL_START_TYPE  == "initial" ) then

  ./xmlchange  -file env_run.xml  -id  RUN_TYPE           -val 'startup'
  ./xmlchange  -file env_run.xml  -id  CONTINUE_RUN       -val 'FALSE'

else if ( $MODEL_START_TYPE  == "continue" ) 

  ./xmlchange  -file env_run.xml  -id  CONTINUE_RUN       -val 'TRUE'

else if ( $MODEL_START_TYPE  == "branch" || $MODEL_START_TYPE  == "hybrid" )

  ./xmlchange  -file env_run.xml  -id  RUN_TYPE           -val $MODEL_START_TYPE
  ./xmlchange  -file env_run.xml  -id  GET_REFCASE        -val $GET_REFCASE
  ./xmlchange  -file env_run.xml  -id  RUN_REFDIR         -val $RUN_REFDIR
  ./xmlchange  -file env_run.xml  -id  RUN_REFCASE        -val $RUN_REFCASE
  ./xmlchange  -file env_run.xml  -id  RUN_REFDATE        -val $RUN_REFDATE
  ./xmlchange  -file env_run.xml  -id  CONTINUE_RUN       -val 'FALSE'
 
else 

  echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
  exit 380

endif 

cat <<EOF >! user_nl_eam
&camexp
!...................
! conditional diag
!...................
metric_name       = 'RHI',      'RHI',
metric_nver       =  72,         72
metric_cmpr_type  =  1,          1
metric_threshold  =  125,        -1
cnd_eval_chkpt    =  'CLDMAC01', 'CLDMAC01'
cnd_end_chkpt     =  'PBCDIAG',  'PBCDIAG'
!
qoi_chkpt   = 'PBCDIAG', 'RAD', 'PACEND','DYNEND','DEEPCU',
              'CLDMAC01','CLDMIC01'
              'CLDMAC02','CLDMIC02'
              'CLDMAC03','CLDMIC03'
              'CLDMAC04','CLDMIC04'
              'CLDMAC05','CLDMIC05'
              'CLDMAC06','CLDMIC06'
!
!
qoi_name = 'RHI', 'Q', 'QSATI'
qoi_nver =  72,    72,  72
!
l_output_state = .true.
l_output_incrm = .true.
!
hist_tape_with_all_output = 1
!
!.......................................................
! history files
!.......................................................
 history_amwg        = .false.
 history_aero_optics = .false.
 history_aerosol     = .false.
!...................
! change init data
!...................
!ncdata = '/compyfs/sunj695/csmruns/init/constance_FC5AV1C-04P2_ne30_ne30_E3SMv1.0_CLIMO_PD_TUNE.cam.i.0000-10-01-00000.nc'
/
EOF

cat <<EOB >! user_nl_elm
&clm_inparm
! finidat = '/compyfs/sunj695/csmruns/init/constance_FC5AV1C-04P2_ne30_ne30_E3SMv1.0_CLIMO_PD_TUNE.clm2.r.0000-10-01-00000.nc'
 check_finidat_fsurdat_consistency = .false.
/
EOB

####Run the model 
./case.submit

endif
