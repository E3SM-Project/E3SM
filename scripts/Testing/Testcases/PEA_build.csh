#!/bin/csh -f
setenv CIMEROOT `./xmlquery CIMEROOT    --value`

./Tools/check_lockedfiles || exit -1

# NOTE - assume that are already in $CASEROOT here
# NOTE - assume that Macros change with differences in MPILIB

# Reset all previous settings
if ( -e env_build.xml.1 )  then
  cp -f env_build.xml.1 env_build.xml
else
  cp -f env_build.xml   env_build.xml.1
endif

set CASE    = `./xmlquery CASE	  --value`
set EXEROOT = `./xmlquery EXEROOT --value`

if(-e mpilib.original) then
  set mpiliborig=`cat mpilib.original`
  ./xmlchange --file env_build.xml --id MPILIB --val $mpiliborig
  ./case.setup -clean
  ./case.setup
else
  set MPILIB  = `./xmlquery MPILIB  --value`
  echo $MPILIB > mpilib.original
endif

#----------------------------
# Build with mpi
#----------------------------

./case.clean_build
./case.build --testmode
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/${CIME_MODEL}.exe $EXEROOT/${CIME_MODEL}.exe.1
cp -f env_build.xml      env_build.xml.1
mv -f Macros Macros.1

#----------------------------
# Build with mpi-serial
#----------------------------

./xmlchange --file env_build.xml --id MPILIB --val mpi-serial

./case.setup -clean 
./case.setup

./case.clean_build
./case.build --testmode
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/${CIME_MODEL}.exe $EXEROOT/${CIME_MODEL}.exe.2
cp -f env_build.xml env_build.xml.2
mv -f Macros Macros.2


