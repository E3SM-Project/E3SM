#!/bin/csh

# script name: setup.testCOSP2.cori.csh

# create and build model

  set PROJECT=acme

# set E3SMROOT to the root to E3SM that contains cime and components subdir

  set E3SMROOT=~/E3SMv1/E3SM_Nov19/E3SM

  echo E3SMroot is $E3SMROOT

  cd $E3SMROOT/cime/scripts

  set compset=FC5AV1C-04P2
  set machine=cori-knl
  set grid=ne30_ne30

# the RUNDIR would be $CSCRATCH/acme_scratch/$CASEID/run

  set CASEID=testCOSPv1_Nov_cosp2

  set CASEDIR=~/ACME/Cases/testing/$machine/$CASEID   #make sure the basepath exist

  ./create_newcase -case $CASEDIR -mach $machine -project $PROJECT -compset $compset -res $grid

  cd $CASEDIR

  ./case.setup > case.setup.log

  ./xmlchange -file env_build.xml CAM_CONFIG_OPTS="-verbose -cosp" -append


  ./xmlchange -file env_run.xml RUN_STARTDATE="0001-01-01",STOP_N="6",STOP_OPTION="nmonths",REST_OPTION="nmonths",REST_N="3",SAVE_TIMING="FALSE",DOUT_S="FALSE",PIO_NUMTASKS="41"

  ./xmlchange JOB_WALLCLOCK_TIME="6:00:00"

# if to build and submit in one script, uncomment the following two lines

# ./case.build
# ./case.submit

  exit
