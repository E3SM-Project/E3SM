#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set xmlquery_data=`./xmlquery -s JGFSEP CASE EXEROOT`
set CASE=`echo $xmlquery_data | awk -F'JGFSEP' '{print $1}'`
set EXEROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $2}'`

./xmlchange -file env_build.xml -id USE_ESMF_LIB   -val TRUE

#------------------------------------------------------------
./xmlchange -file env_build.xml -id COMP_INTERFACE -val MCT
./$CASE.clean_build

./$CASE.build
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.mct
cp -f env_build.xml      env_build.xml.mct

#------------------------------------------------------------
./xmlchange -file env_build.xml -id COMP_INTERFACE -val ESMF
./$CASE.clean_build

./$CASE.build
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.esmf
cp -f env_build.xml  env_build.xml.esmf




