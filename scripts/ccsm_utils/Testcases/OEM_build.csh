#!/bin/csh -f 
./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cd $CASEROOT

./cesm_setup -clean
./cesm_setup 
./xmlchange -file env_run.xml -id BFBFLAG -val TRUE
echo "b4b_flag=.true." >> user_nl_pop2
./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.oem1
cp -f env_build.xml env_build.xml.oem1
cp -f env_mach_pes.xml env_mach_pes.xml.oem1

source ./Tools/ccsm_getenv || exit -3

if ( $NTASKS_ATM > 1) then
	@ ntask = $NTASKS_ATM / 2
	./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $ntask
endif
if ($NTHRDS_ATM == 1) then
	@ nthrd = $NTHRDS_ATM * 2
	./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthrd
endif
if ( $NTASKS_LND > 1) then
	@ ntask = $NTASKS_LND / 2
	./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $ntask
endif
if ($NTHRDS_LND == 1) then
	@ nthrd = $NTHRDS_LND * 2
	./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthrd
endif
if ( $NTASKS_OCN > 1) then
	@ ntask = $NTASKS_OCN / 2
	./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $ntask
endif
if ($NTHRDS_OCN == 1) then
	@ nthrd = $NTHRDS_OCN * 2
	./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $nthrd
endif
if ( $NTASKS_ICE > 1) then
	@ ntask = $NTASKS_ICE / 2
	./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $ntask
endif
if ($NTHRDS_ICE == 1) then
	@ nthrd = $NTHRDS_ICE * 2
	./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthrd
endif
if ( $NTASKS_ICE > 1) then
	@ ntask = $NTASKS_ICE / 2
	./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $ntask
endif
if ($NTHRDS_ICE == 1) then
	@ nthrd = $NTHRDS_ICE * 2
	./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthrd
endif
if ( $NTASKS_GLC > 1) then
	@ ntask = $NTASKS_GLC / 2
	./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $ntask
endif
if ($NTHRDS_ICE == 1) then
	@ nthrd = $NTHRDS_GLC * 2
	./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $nthrd
endif

./cesm_setup -clean
./cesm_setup
./$CASE.clean_build -all 
./$CASE.build

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.oem2
cp -f env_build.xml env_build.xml.oem2
cp -f env_mach_pes.xml env_mach_pes.xml.oem2

cp -f env_build.xml.oem1 env_build.xml
cp -f env_mach_pes.xml.oem1 env_mach_pes.xml
./cesm_setup -clean
./cesm_setup
cp -f env_build.xml.oem1 env_build.xml
cp -f env_mach_pes.xml.oem1 env_mach_pes.xml
cp -f env_mach_pes.xml LockedFiles/env_mach_pes.xml.locked
cp -f env_build.xml LockedFiles/env_build.xml.locked
