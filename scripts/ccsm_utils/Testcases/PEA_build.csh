#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cd $CASEROOT

./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.mpi
cp -f env_build.xml env_build.xml.mpi

./xmlchange -file env_build.xml -id MPILIB -val mpi-serial
./$CASE.clean_build
./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.mpi-serial
cp -f env_build.xml env_build.xml.mpi-serial

cp -f env_build.xml.mpi env_build.xml
cp -f env_build.xml LockedFiles/env_build.xml.locked

