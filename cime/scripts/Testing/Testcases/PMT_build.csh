#!/bin/csh -f 

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE        = `./xmlquery CASE		-value`
set EXEROOT     = `./xmlquery EXEROOT		-value`

# Reset all previous env_mach_pes settings
if ( -e env_mach_pes.xml.1 )  then
  cp -f env_mach_pes.xml.1 env_mach_pes.xml
endif

./cesm_setup -clean -testmode
./cesm_setup 

cp -f env_mach_pes.xml env_mach_pes.xml.1
cp -f env_mach_pes.xml LockedFiles/env_mach_pes.xml.locked

./xmlchange -file env_run.xml -id BFBFLAG -val TRUE
echo "b4b_flag=.true." >> user_nl_pop

./$CASE.build
if ($status != 0) then
   echo "Error: build for default PE layout failed" >! ./TestStatus
   echo "CFAIL $CASE" > ./TestStatus
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.1
cp -f env_build.xml    env_build.xml.1

set NTASKS_ATM  = `./xmlquery NTASKS_ATM	-value`
set NTASKS_LND  = `./xmlquery NTASKS_LND	-value`
set NTASKS_ROF  = `./xmlquery NTASKS_ROF	-value`
set NTASKS_WAV  = `./xmlquery NTASKS_WAV	-value`
set NTASKS_OCN  = `./xmlquery NTASKS_OCN	-value`
set NTASKS_ICE  = `./xmlquery NTASKS_ICE	-value`
set NTASKS_GLC  = `./xmlquery NTASKS_GLC	-value`
set NTASKS_CPL  = `./xmlquery NTASKS_CPL	-value`

set NTHRDS_ATM  = `./xmlquery NTHRDS_ATM	-value`
set NTHRDS_LND  = `./xmlquery NTHRDS_LND	-value`
set NTHRDS_ROF  = `./xmlquery NTHRDS_ROF	-value`
set NTHRDS_WAV  = `./xmlquery NTHRDS_WAV	-value`
set NTHRDS_OCN  = `./xmlquery NTHRDS_OCN	-value`
set NTHRDS_ICE  = `./xmlquery NTHRDS_ICE	-value`
set NTHRDS_GLC  = `./xmlquery NTHRDS_GLC	-value`
set NTHRDS_CPL  = `./xmlquery NTHRDS_CPL	-value`

# Halve the number of tasks 
if ( $NTASKS_ATM > 1) then
  @ ntask = $NTASKS_ATM / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $ntask
endif
if ( $NTASKS_LND > 1) then
    @ ntask = $NTASKS_LND / 2
    ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $ntask
endif
if ( $NTASKS_ICE > 1) then
    @ ntask = $NTASKS_ICE / 2
    ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $ntask
endif
if ( $NTASKS_OCN > 1) then
    @ ntask = $NTASKS_OCN / 2
    ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $ntask
endif
if ( $NTASKS_GLC > 1) then
    @ ntask = $NTASKS_GLC / 2
    ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $ntask
endif

# Double the number of threads
if ($NTHRDS_ATM == 1) then
    @ nthrd = $NTHRDS_ATM * 2
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthrd
endif
if ($NTHRDS_LND == 1) then
    @ nthrd = $NTHRDS_LND * 2
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthrd
endif
if ($NTHRDS_ICE == 1) then
    @ nthrd = $NTHRDS_ICE * 2
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthrd
endif
if ($NTHRDS_OCN == 1) then
    @ nthrd = $NTHRDS_OCN * 2
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $nthrd
endif
if ($NTHRDS_ICE == 1) then
    @ nthrd = $NTHRDS_ICE * 2
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthrd
endif
if ($NTHRDS_ICE == 1) then
    @ nthrd = $NTHRDS_GLC * 2
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $nthrd
endif

# Build with half the tasks and double the threads
./cesm_setup -clean -testmode
./cesm_setup
./xmlchange -file env_build.xml -id SMP_BUILD -val 0

./cesm_setup -clean
./cesm_setup
./$CASE.clean_build -all 

./$CASE.build
if ($status != 0) then
   echo "Error: build for default half tasks twice threads failed" >! ./TestStatus
   echo "CFAIL $CASE" > ./TestStatus
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.2
cp -f env_build.xml env_build.xml.2
cp -f env_mach_pes.xml env_mach_pes.xml.2
