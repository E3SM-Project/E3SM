#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cp -f env_mach_pes.xml env_mach_pes.xml.nck1
cp -f env_build.xml env_build.xml.nck1

cd $CASEROOT
./$CASE.clean_build

./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_WAV  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_OCN  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_ICE  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_GLC  -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM  -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND  -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF  -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV  -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN  -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE  -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC  -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL  -val 0
if ( $NTASKS_ATM > 1 ) then
  @ ntask = $NTASKS_ATM / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM  -val $ntask
endif
if ( $NTASKS_LND > 1 ) then
  @ ntask = $NTASKS_LND / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_LND  -val $ntask
endif
if ( $NTASKS_ROF > 1 ) then
  @ ntask = $NTASKS_ROF / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF  -val $ntask
endif
if ( $NTASKS_WAV > 1 ) then
  @ ntask = $NTASKS_WAV / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV  -val $ntask
endif
if ( $NTASKS_OCN > 1 ) then
  @ ntask = $NTASKS_OCN / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN  -val $ntask
endif
if ( $NTASKS_ICE > 1 ) then
  @ ntask = $NTASKS_ICE / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE  -val $ntask
endif
if ( $NTASKS_GLC > 1 ) then
  @ ntask = $NTASKS_GLC / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC  -val $ntask
endif
./xmlchange -file env_build.xml -id NINST_BUILD -val 0
./cesm_setup -clean
./cesm_setup
./$CASE.build

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.nck1  || exit -9
cp -f env_mach_pes.xml env_mach_pes.xml.nck1
cp -f env_build.xml env_build.xml.nck1

./$CASE.clean_build 

source ./Tools/ccsm_getenv     || exit -2
./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_WAV  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_OCN  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_ICE  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_GLC  -val 2
@ rootp = 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL  -val $rootp
@ ntask = $NTASKS_ATM * 2
@ rootp = $rootp + $ntask
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM  -val $ntask
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM  -val $rootp
@ ntask = $NTASKS_LND * 2
@ rootp = $rootp + $ntask
./xmlchange -file env_mach_pes.xml -id NTASKS_LND  -val $ntask
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND  -val $rootp
@ ntask = $NTASKS_ROF * 2
@ rootp = $rootp + $ntask
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF  -val $ntask
./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF  -val $rootp
@ ntask = $NTASKS_WAV * 2
@ rootp = $rootp + $ntask
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV  -val $ntask
./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV  -val $rootp
@ ntask = $NTASKS_OCN * 2
@ rootp = $rootp + $ntask
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN  -val $ntask
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN  -val $rootp
@ ntask = $NTASKS_ICE * 2
@ rootp = $rootp + $ntask
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE  -val $ntask
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE  -val $rootp
@ ntask = $NTASKS_GLC * 2
@ rootp = $rootp + $ntask
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC  -val $ntask
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC  -val $rootp

./xmlchange -file env_build.xml -id NINST_BUILD -val 0
./cesm_setup -clean
./cesm_setup
./$CASE.build

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.nck2  || exit -9
cp -f env_mach_pes.xml env_mach_pes.xml.nck2
cp -f env_build.xml env_build.xml.nck2

cp -f env_mach_pes.xml.nck1 env_mach_pes.xml
cp -f env_build.xml.nck1 env_build.xml
cp -f env_mach_pes.xml LockedFiles/env_mach_pes.xml.locked
cp -f env_build.xml LockedFiles/env_build.xml.locked

