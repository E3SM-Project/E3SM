#!/bin/csh -f

####################################################
# Set case variables
####################################################
setenv cesmsrc /global/u1/m/mickelso/cesm1_3_alpha04b/
setenv res ne30_g16
setenv compset B1850C5
setenv mach edison_intel
setenv casedir $SCRATCH/tests/cesm_timing_tests/B1850C5.ne30_g16/
setenv casestr  _B1850C5_ne30_g16__
setenv run_len 10

# Select either FV or SE below
#setenv DYCORE FV
setenv DYCORE SE

####################################################
# Set the location of the load balancing results
####################################################
setenv results_dir /global/u1/m/mickelso/neos-version-example/results/

####################################################
# Set the test layouts to produce the scaling curves
####################################################
setenv NTHRDS_VAL 1

# Set the Task Counts
setenv TASK_ATM  "128,256,512,1024,2048"
setenv TASK_LND  "128,256,512,1024,2048"
setenv TASK_ROF  "128,256,512,1024,2048"
setenv TASK_ICE  "128,256,512,1024,2048"
setenv TASK_OCN  "128,256,512,1024,2048"
setenv TASK_CPL  "128,256,512,1024,2048"
setenv TASK_WAV  "128,256,512,1024,2048"
setenv TASK_GLC  "1,1,1,1,1"

# Set Root Locations
setenv ROOT_ATM  "0,0,0,0,0"
setenv ROOT_LND  "0,0,0,0,0"
setenv ROOT_ROF  "0,0,0,0,0"
setenv ROOT_ICE  "0,0,0,0,0"
setenv ROOT_OCN  "0,0,0,0,0"
setenv ROOT_CPL  "0,0,0,0,0"
setenv ROOT_WAV  "0,0,0,0,0"
setenv ROOT_GLC  "0,0,0,0,0"

####################################################
# Set the target task counts (ATM(LND+ICE) + OCN)
####################################################
setenv TARGET_TASKS "256,512,1024"

