#!/bin/csh

 

  set E3SMROOT=/qfs/people/shpu881/E3SM/E3SM_alpha5_59_v2candidate_NGD_Conv

  #

  set compset=F2010SC5-P3   # Or you may use F20TR, or F2010

  set machine=compy

  set grid=ne30pg2_r05_oECv3

  set compiler=intel

 

  set PROJECT=e3sm

 

 

  set CASEID=2021April14_Test.$compset

 

  set SCRATCH=/compyfs/shpu881/E3SM_simulations/alpha5_59_v2candidate_NGD_Conv

  set CASEDIR=$SCRATCH/E3SM_testings/$CASEID

  set CASE_SCRIPTS=$CASEDIR/case_scripts

  set EXEROOT=$CASEDIR/bld

  set RUNDIR=$CASEDIR/run

 

# create a new case

  rm -rf $CASE_SCRIPTS
  rm -rf $EXEROOT
  rm -rf $RUNDIR

  cd $E3SMROOT/cime/scripts

  ./create_newcase -case $CASEID -mach $machine -project $PROJECT -compset $compset -res $grid --script-root $CASE_SCRIPTS


  cd $CASE_SCRIPTS

 

  ./xmlchange EXEROOT=$EXEROOT

  ./xmlchange RUNDIR=$RUNDIR

  ./xmlchange CAM_TARGET=theta-l  # for theta-l + SL tracer transport

  ./xmlchange --id CAM_CONFIG_OPTS --val "-mach $machine -phys default -clubb_sgs -microphys p3 -chem linoz_mam4_resus_mom_soag -rain_evap_to_coarse_aero -nlev 72 -bc_dep_to_snow_updates -cosp -cosp"

# turn cosp on if plan to use the simulation for other purpose as well
#
  ./xmlchange CAM_CONFIG_OPTS="-cosp" --append

# Override default pe-layout
#
#  cp /global/cfs/cdirs/e3sm/wlin/share/FY20-Q4-MCS/env_mach_pes.xml.169-nodes env_mach_pes.xml
   cp /qfs/people/shpu881/E3SM/E3SM_alpha5_59_v2candidate_NGD_Conv/env_mach_pes_files/1200x1.env_mach_pes.xml env_mach_pes.xml  

  cat >> user_nl_eam <<EOF

  cosp_lite = .true.

  avgflag_pertape = 'A','A','A','I','A'

  nhtfrq=-24,-24,-6,-3

  fincl2 = 'PRECT','U200','V200'

  fincl3 = 'OMEGA500','PRECT','U200','U850','FLUT','UBOT:I','VBOT:I','Z700:I','PSL:I','PS:I','TREFHT','TGCLDLWP','TGCLDIWP'

  fincl4='PRECT:A','LHFLX:A','SHFLX:A'

  mfilt=1,30,120,240,720


 state_debug_checks = .true.
 history_budget = .true.
 

EOF

 

  ./case.setup --reset > case.setup.log

 
  ./xmlchange --id DEBUG --val 'true'
  ./xmlchange --id JOB_QUEUE --val 'short' #'debug'
  ./xmlchange JOB_WALLCLOCK_TIME="01:00:00" 
  ./xmlchange -file env_run.xml RUN_STARTDATE="2010-01-01",STOP_N="24",STOP_OPTION="nmonths",REST_OPTION="nmonths",REST_N="3",SAVE_TIMING="TRUE",DOUT_S="FALSE"

  

  ./case.build >& case.log-build

  ./case.submit
