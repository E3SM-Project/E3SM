#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

# NOTE - assume that are already in $CASEROOT here
# NOTE - assume that Macros change with differences in MPILIB

# Reset all previous settings
if ( -e env_build.xml.1 )  then
  cp -f env_build.xml.1 env_build.xml
else
  cp -f env_build.xml   env_build.xml.1
endif

set xmlquery_data=`./xmlquery -s JGFSEP CASE EXEROOT MPILIB`
set CASE=`echo $xmlquery_data | awk -F'JGFSEP' '{print $1}'`
set EXEROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $2}'`
set MPILIB=`echo $xmlquery_data | awk -F'JGFSEP' '{print $3}'`

if(-e mpilib.original) then
  set mpiliborig=`cat mpilib.original`
  ./xmlchange -file env_build.xml -id MPILIB -val $mpiliborig
  ./cesm_setup -clean
  ./cesm_setup
else
  echo $MPILIB > mpilib.original
endif

#----------------------------
# Build with mpi
#----------------------------

./$CASE.clean_build
./$CASE.build
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.1
cp -f env_build.xml      env_build.xml.1
mv -f Macros Macros.1

#----------------------------
# Build with mpi-serial
#----------------------------

./xmlchange -file env_build.xml -id MPILIB -val mpi-serial

./cesm_setup -clean 
./cesm_setup

./$CASE.clean_build
./$CASE.build
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.2
cp -f env_build.xml env_build.xml.2
mv -f Macros Macros.2


