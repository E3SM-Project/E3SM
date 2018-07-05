#!/bin/bash

#######################################################################
#######################################################################
#######  Script to run ACME in ensemble mode using CIME's multi instance feature
######   Currently runs a test FC5 with 2 active cases
#######################################################
#######  BEGIN USER DEFINED SETTINGS

# Set the name of your case here
casename=ICLM45GLCMEC

# Set the case directory here
casedirectory=$PROJWORK/cli115/$USER
 
# Name of machine you are running on (i.e. edison, anvil, etc)                                                    
machine=titan

# Name of project to run on, if submitting to queue
projectname=cli115

# Name of run type
compset=ICLM45GLCMEC

# grid resolution
grid=ne30_ne30

# User enter any needed modules to load or use below
module load python/2.7.9

# Directory where code lives
code_dir=$HOME/ACME

# Code tag name 
code_tag=ACME   
                                                         
# Create new case
$code_dir/$code_tag/cime/scripts/create_newcase -case $casedirectory/$casename -mach $machine -project $projectname -compset $compset -res $grid 

cd $casedirectory/$casename

./pelayout 
 
./xmlquery TOTALPES

# Enter CAM namelist options 
#cat <<EOF >> user_nl_cam
# fincl1 = 'SWCF:A','LWCF:A'
# hist_nhtfrq = 0,-24
#EOF
# Enter CLM namelist options 
cat <<EOF >> user_nl_clm
 flndtopo = '/lustre/atlas1/cli900/world-shared/cesm/inputdata/lnd/clm2/griddata/topodata_0.9x1.25_USGS_070110.nc'
 fglcmask = '/lustre/atlas1/cli900/world-shared/cesm/inputdata/lnd/clm2/griddata/glcmaskdata_0.9x1.25_GIS_AIS.nc'
EOF
cat <<EOF >> user_nl_mali
 mali_use_albany='FALSE'
EOF

# set up the case
  ./case.setup

# Build the case 
  ./case.build

  ./xmlchange STOP_N=32
  ./xmlchange JOB_WALLCLOCK_TIME=02:30:00

# Submit case 
  ./case.submit
  
