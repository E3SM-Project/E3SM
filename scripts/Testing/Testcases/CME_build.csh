#!/bin/csh -f
setenv CIMEROOT `./xmlquery CIMEROOT    -value`

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE     = `./xmlquery CASE    -value`
set EXEROOT  = `./xmlquery EXEROOT -value`

./xmlchange USE_ESMF_LIB=TRUE

#------------------------------------------------------------
./xmlchange COMP_INTERFACE=MCT
./case.clean_build

./case.build
if ($status != 0) then
   echo "Error: build for MCT failed" >! ./TestStatus
   echo "CFAIL $CASE" > ./TestStatus
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.mct
cp -f env_build.xml      env_build.xml.mct

#------------------------------------------------------------
./xmlchange COMP_INTERFACE=ESMF
./case.clean_build

./case.build
if ($status != 0) then
   echo "Error: build for ESMF failed" >! ./TestStatus
   echo "CFAIL $CASE" > ./TestStatus
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.esmf
cp -f env_build.xml  env_build.xml.esmf




