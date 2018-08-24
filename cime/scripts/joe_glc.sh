#!/usr/bin/env bash


# -----
# TITAN
# -----
# module load python
#
# export CHARGE_ACCOUNT=cli106ice
#
# case=$PROJWORK/cli106/$USER/snow_joe_0
#
# ./create_newcase -case $case -compset IGCLM45_MLI -res f09_g16_a -project cli106
#
# cd $case
#
# ./pelayout
#
# ./xmlquery TOTALPES
#
# ./case.setup
#
# ./case.build
#
# ./xmlchange STOP_N=32
# ./xmlchange JOB_WALLCLOCK_TIME=02:30:00
#
# ./case.submit


# ------
# EDISON
# ------
case=$CSCRATCH/kennedyj/snow_joe_0

./create_newcase -case $case -compset IGCLM45_MLI -res f09_g16_a -project m1795

cd $case

./pelayout

./xmlquery TOTALPES

./case.setup

./case.build

./xmlchange STOP_N=32
./xmlchange JOB_WALLCLOCK_TIME=02:30:00

./case.submit
