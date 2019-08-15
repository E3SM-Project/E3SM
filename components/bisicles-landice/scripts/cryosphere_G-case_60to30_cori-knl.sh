#####################################################################
# Template for E3SM Cryosphere science campaign simulations: G-cases.
# Written for NERSC, directory modifications needed for other LCFs.
# Author: Darin Comeau, LANL.
#
# Things to check:
# - change USER
# - change REPO_ROOT, CASE_ROOT directories if desired
# - change REPO_NAME if an existing repository is desired
# - check COMPSET
# - check GRID
# - check add a tag to E3SM_CASE if desired
# - check PE_LAYOUT
# - check RUN_OPTIONS, default is 5 day test
# 
#####################################################################

#####################################################################
# Define user paths
#####################################################################

export USER=dcomeau
export PROJECT=m2833
export MACHINE=cori-knl
export REPO_ROOT=/project/projectdirs/$PROJECT/$USER/e3sm_repos
export CASE_ROOT=/project/projectdirs/$PROJECT/$USER/e3sm_cases
export RUN_ROOT=/global/cscratch1/sd/$USER/acme_scratch/$MACHINE

if [ ! -d "$REPO_ROOT" ]; then
	mkdir -p $REPO_ROOT
fi
if [ ! -d "$CASE_ROOT" ]; then
	mkdir -p $CASE_ROOT
fi

#####################################################################
# Download repository
#####################################################################

# export REPO_NAME=maint-1.2.`date +"%Y%m%d"`
export REPO_NAME=maint-1.2
export E3SM_ROOT=$REPO_ROOT/$REPO_NAME/E3SM

if [ ! -d "$REPO_ROOT/$REPO_NAME" ]; then
	mkdir -p $REPO_ROOT/$REPO_NAME
	cd $REPO_ROOT/$REPO_NAME
	git clone git@github.com:E3SM-Project/E3SM.git
	cd E3SM
	git checkout origin/maint-1.2
	git submodule update --init
else
	cd $REPO_ROOT/$REPO_NAME/E3SM
	git fetch
	git reset --hard origin/maint-1.2
	git submodule update --init	
fi

#####################################################################
# Define case:
# COMPSET options:
# GMPAS-IAF, GMPAS-IAF-ISMF, GMPAS-DIB-IAF, GMPAS-DIB-IAF-ISMF
# GRID options:
# T62_oEC60to30v3wLI, T62_oRRS30to10v3wLI
#####################################################################

export COMPSET=GMPAS-DIB-IAF-ISMF
export GRID=T62_oEC60to30v3wLI
export COMPILER=intel
export E3SM_CASE=`date +"%Y%m%d"`.${COMPSET}.${GRID}.${MACHINE}

#####################################################################
# Create case
#####################################################################

cd $E3SM_ROOT/cime/scripts
./create_newcase \
-case $CASE_ROOT/$E3SM_CASE \
-compiler $COMPILER \
-mach $MACHINE \
-project $PROJECT \
-compset $COMPSET \
-res $GRID

cd $CASE_ROOT/$E3SM_CASE

#####################################################################
# PE_LAYOUT
#####################################################################

## Medium layout
# export NPROCS_CPL=640
# export NPROCS_ATM=640
# export NPROCS_LND=640
# export NPROCS_ROF=640
# export NPROCS_ICE=640
# export NPROCS_OCN=960
# export NPROCS_GLC=640
# export NPROCS_WAV=640

# export NTHRDS_CPL=1
# export NTHRDS_ATM=1
# export NTHRDS_LND=1
# export NTHRDS_ROF=1
# export NTHRDS_ICE=1
# export NTHRDS_OCN=1
# export NTHRDS_GLC=1
# export NTHRDS_WAV=1

# export ROOTPE_CPL=0
# export ROOTPE_ATM=0
# export ROOTPE_LND=0
# export ROOTPE_ROF=0
# export ROOTPE_ICE=0
# export ROOTPE_OCN=680
# export ROOTPE_GLC=0
# export ROOTPE_WAV=0

## Large layout
# export NPROCS_CPL=1200
# export NPROCS_ATM=1200
# export NPROCS_LND=1200
# export NPROCS_ROF=1200
# export NPROCS_ICE=1200
# export NPROCS_OCN=1600
# export NPROCS_GLC=1200
# export NPROCS_WAV=1200

# export NTHRDS_CPL=1
# export NTHRDS_ATM=1
# export NTHRDS_LND=1
# export NTHRDS_ROF=1
# export NTHRDS_ICE=1
# export NTHRDS_OCN=1
# export NTHRDS_GLC=1
# export NTHRDS_WAV=1

# export ROOTPE_CPL=0
# export ROOTPE_ATM=0
# export ROOTPE_LND=0
# export ROOTPE_ROF=0
# export ROOTPE_ICE=0
# export ROOTPE_OCN=1224
# export ROOTPE_GLC=0
# export ROOTPE_WAV=0

#####################################################################

export NPROCS_CPL=1200
export NPROCS_ATM=1200
export NPROCS_LND=1200
export NPROCS_ROF=1200
export NPROCS_ICE=1200
export NPROCS_OCN=1600
export NPROCS_GLC=1200
export NPROCS_WAV=1200

export NTHRDS_CPL=1
export NTHRDS_ATM=1
export NTHRDS_LND=1
export NTHRDS_ROF=1
export NTHRDS_ICE=1
export NTHRDS_OCN=1
export NTHRDS_GLC=1
export NTHRDS_WAV=1

export ROOTPE_CPL=0
export ROOTPE_ATM=0
export ROOTPE_LND=0
export ROOTPE_ROF=0
export ROOTPE_ICE=0
export ROOTPE_OCN=1224
export ROOTPE_GLC=0
export ROOTPE_WAV=0

./xmlchange -file env_mach_pes.xml  -id NTASKS_CPL  -val $NPROCS_CPL
./xmlchange -file env_mach_pes.xml  -id NTASKS_ATM  -val $NPROCS_ATM
./xmlchange -file env_mach_pes.xml  -id NTASKS_LND  -val $NPROCS_LND
./xmlchange -file env_mach_pes.xml  -id NTASKS_ROF  -val $NPROCS_ROF
./xmlchange -file env_mach_pes.xml  -id NTASKS_ICE  -val $NPROCS_ICE
./xmlchange -file env_mach_pes.xml  -id NTASKS_OCN  -val $NPROCS_OCN
./xmlchange -file env_mach_pes.xml  -id NTASKS_GLC  -val $NPROCS_GLC
./xmlchange -file env_mach_pes.xml  -id NTASKS_WAV  -val $NPROCS_WAV

./xmlchange -file env_mach_pes.xml  -id NTHRDS_CPL  -val $NTHRDS_CPL
./xmlchange -file env_mach_pes.xml  -id NTHRDS_ATM  -val $NTHRDS_ATM
./xmlchange -file env_mach_pes.xml  -id NTHRDS_LND  -val $NTHRDS_LND
./xmlchange -file env_mach_pes.xml  -id NTHRDS_ROF  -val $NTHRDS_ROF
./xmlchange -file env_mach_pes.xml  -id NTHRDS_ICE  -val $NTHRDS_ICE
./xmlchange -file env_mach_pes.xml  -id NTHRDS_OCN  -val $NTHRDS_OCN
./xmlchange -file env_mach_pes.xml  -id NTHRDS_GLC  -val $NTHRDS_GLC
./xmlchange -file env_mach_pes.xml  -id NTHRDS_WAV  -val $NTHRDS_WAV

./xmlchange -file env_mach_pes.xml  -id ROOTPE_CPL  -val $ROOTPE_CPL
./xmlchange -file env_mach_pes.xml  -id ROOTPE_ATM  -val $ROOTPE_ATM
./xmlchange -file env_mach_pes.xml  -id ROOTPE_LND  -val $ROOTPE_LND
./xmlchange -file env_mach_pes.xml  -id ROOTPE_ROF  -val $ROOTPE_ROF
./xmlchange -file env_mach_pes.xml  -id ROOTPE_ICE  -val $ROOTPE_ICE
./xmlchange -file env_mach_pes.xml  -id ROOTPE_OCN  -val $ROOTPE_OCN
./xmlchange -file env_mach_pes.xml  -id ROOTPE_GLC  -val $ROOTPE_GLC
./xmlchange -file env_mach_pes.xml  -id ROOTPE_WAV  -val $ROOTPE_WAV

#####################################################################
# Setup and build
#####################################################################

./case.setup
./case.build

#####################################################################
# RUN_OPTIONS
#
# SMOKE TEST:
# ./xmlchange -file env_run.xml -id STOP_OPTION -val ndays
# ./xmlchange -file env_run.xml -id STOP_N -val 5
# ./xmlchange -file env_run.xml -id REST_OPTION -val nmonths
# ./xmlchange -file env_run.xml -id REST_N -val 1
# ./xmlchange -file env_batch.xml -subgroup case.run -id JOB_QUEUE -val debug
# ./xmlchange -file env_batch.xml -subgroup case.run -id JOB_WALLCLOCK_TIME -val "00:30:00"
#
# FULL RUN:
# ./xmlchange -file env_run.xml -id STOP_OPTION -val nyears
# ./xmlchange -file env_run.xml -id STOP_N -val 3
# ./xmlchange -file env_run.xml -id REST_OPTION -val nyears
# ./xmlchange -file env_run.xml -id REST_N -val 1
# ./xmlchange -file env_batch.xml -subgroup case.run -id JOB_QUEUE -val regular
# ./xmlchange -file env_batch.xml -subgroup case.run -id JOB_WALLCLOCK_TIME -val "12:00:00"
#
# AFTER FIRST FULL RUN:
# ./xmlchange -file env_run.xml -id CONTINUE_RUN -val TRUE
# ./xmlchange -file env_run.xml -id RESUBMIT -val 12
#
#####################################################################

./xmlchange -file env_run.xml -id STOP_OPTION -val ndays
./xmlchange -file env_run.xml -id STOP_N -val 5
./xmlchange -file env_run.xml -id REST_OPTION -val nmonths
./xmlchange -file env_run.xml -id REST_N -val 1
./xmlchange -file env_batch.xml -subgroup case.run -id JOB_QUEUE -val debug
./xmlchange -file env_batch.xml -subgroup case.run -id JOB_WALLCLOCK_TIME -val "00:30:00"

./case.submit

#####################################################################
# Open up permissions
#####################################################################

cd $RUN_ROOT
chgrp acme ${E3SM_CASE}
chmod 750 ${E3SM_CASE}
chmod g+s ${E3SM_CASE}
setfacl -d -m g::rx ${E3SM_CASE}
chgrp -R acme ${E3SM_CASE}
chmod -R g-w ${E3SM_CASE}

#####################################################################
# Move to run directory to check on job
#
# cd $RUN_ROOT/$E3SM_CASE/run
#
#####################################################################
