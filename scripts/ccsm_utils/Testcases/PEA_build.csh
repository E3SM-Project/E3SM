#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cd $CASEROOT
./$CASE.clean_build

if(-e mpilib.original) then
  set mpiliborig=`cat mpilib.original`
  ./xmlchange -file env_build.xml -id MPILIB -val $mpiliborig
  ./cesm_setup -clean
  ./cesm_setup
else
  echo $MPILIB > mpilib.original
endif

./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.mpi
cp -f env_build.xml env_build.xml.mpi
mv -f Macros Macros.mpi

./xmlchange -file env_build.xml -id MPILIB -val mpi-serial
./cesm_setup -clean
./cesm_setup
./$CASE.clean_build
./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.mpi-serial
cp -f env_build.xml env_build.xml.mpi-serial
cp -f Macros Macros.mpi-serial

cp -f env_build.xml.mpi env_build.xml
cp -f env_build.xml LockedFiles/env_build.xml.locked

