#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cd $CASEROOT

./xmlchange -file env_build.xml -id COMP_INTERFACE -val MCT
./xmlchange -file env_build.xml -id USE_ESMF_LIB -val TRUE
./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.mct
cp -f env_build.xml env_build.xml.mct

./xmlchange -file env_build.xml -id COMP_INTERFACE -val ESMF
./$CASE.clean_build
./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.esmf
cp -f env_build.xml env_build.xml.esmf

cp -f env_build.xml.mct env_build.xml
cp -f env_build.xml LockedFiles/env_build.xml.locked


