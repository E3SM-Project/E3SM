#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE        = `./xmlquery CASE     -value`
set EXEROOT     = `./xmlquery EXEROOT  -value`

#------------------------------------------------------------
./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_WAV  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_OCN  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_ICE  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_GLC  -val 1

set NTASKS_ATM  = `./xmlquery NTASKS_ATM -value`
set NTASKS_LND  = `./xmlquery NTASKS_LND -value`
set NTASKS_ROF  = `./xmlquery NTASKS_ROF -value`
set NTASKS_WAV  = `./xmlquery NTASKS_WAV -value`
set NTASKS_OCN  = `./xmlquery NTASKS_OCN -value`
set NTASKS_ICE  = `./xmlquery NTASKS_ICE -value`
set NTASKS_GLC  = `./xmlquery NTASKS_GLC -value`
set NTASKS_CPL  = `./xmlquery NTASKS_CPL -value`

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

./cesm_setup -clean -testmode
./cesm_setup

./$CASE.clean_build

./$CASE.build
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.1  || exit -9
cp -f env_mach_pes.xml   env_mach_pes.xml.1
cp -f env_build.xml      env_build.xml.1

#------------------------------------------------------------
./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_WAV  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_OCN  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_ICE  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_GLC  -val 2

set NTASKS_ATM  = `./xmlquery NTASKS_ATM -value`
set NTASKS_LND  = `./xmlquery NTASKS_LND -value`
set NTASKS_ROF  = `./xmlquery NTASKS_ROF -value`
set NTASKS_WAV  = `./xmlquery NTASKS_WAV -value`
set NTASKS_OCN  = `./xmlquery NTASKS_OCN -value`
set NTASKS_ICE  = `./xmlquery NTASKS_ICE -value`
set NTASKS_GLC  = `./xmlquery NTASKS_GLC -value`
set NTASKS_CPL  = `./xmlquery NTASKS_CPL -value`

@ ntask = $NTASKS_ATM * 2
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM  -val $ntask
@ ntask = $NTASKS_LND * 2
./xmlchange -file env_mach_pes.xml -id NTASKS_LND  -val $ntask
@ ntask = $NTASKS_ROF * 2
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF  -val $ntask
@ ntask = $NTASKS_WAV * 2
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV  -val $ntask
@ ntask = $NTASKS_OCN * 2
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN  -val $ntask
@ ntask = $NTASKS_ICE * 2
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE  -val $ntask
@ ntask = $NTASKS_GLC * 2
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC  -val $ntask

./xmlchange -file env_build.xml -id NINST_BUILD -val 0

./cesm_setup -clean -testmode
./cesm_setup

./$CASE.clean_build 

./$CASE.build
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.2  || exit -9
cp -f env_mach_pes.xml   env_mach_pes.xml.2
cp -f env_build.xml      env_build.xml.2


