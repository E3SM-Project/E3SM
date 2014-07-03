#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cd $CASEROOT

# Build with default PE layout
./cesm_setup -clean
./cesm_setup
./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.prs1
cp -f env_build.xml env_build.xml.prs1
cp -f env_mach_pes.xml env_mach_pes.xml.prs1

# Halve the number of tasks and threads
source ./Tools/ccsm_getenv || exit -3

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
if ( $NTASKS_CPL > 1 ) then
  @ ntask = $NTASKS_CPL / 2
  ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL  -val $ntask
endif
if ( $NTHRDS_ATM > 1 ) then
  @ ntask = $NTHRDS_ATM / 2
  ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM  -val $ntask
endif
if ( $NTHRDS_LND > 1 ) then
  @ ntask = $NTHRDS_LND / 2
  ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND  -val $ntask
endif
if ( $NTHRDS_ROF > 1 ) then
  @ ntask = $NTHRDS_ROF / 2
  ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF  -val $ntask
endif
if ( $NTHRDS_WAV > 1 ) then
  @ ntask = $NTHRDS_WAV / 2
  ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV  -val $ntask
endif
if ( $NTHRDS_OCN > 1 ) then
  @ ntask = $NTHRDS_OCN / 2
  ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN  -val $ntask
endif
if ( $NTHRDS_ICE > 1 ) then
  @ ntask = $NTHRDS_ICE / 2
  ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE  -val $ntask
endif
if ( $NTHRDS_GLC > 1 ) then
  @ ntask = $NTHRDS_GLC / 2
  ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC  -val $ntask
endif
if ( $NTHRDS_CPL > 1 ) then
  @ ntask = $NTHRDS_CPL / 2
  ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL  -val $ntask
endif

# Build with half the tasks and threads
./cesm_setup -clean
./cesm_setup
./$CASE.build
mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.prs2
cp -f env_build.xml env_build.xml.prs2
cp -f env_mach_pes.xml env_mach_pes.xml.prs2

cp -f env_build.xml.prs1 env_build.xml
cp -f env_mach_pes.xml.prs1 env_mach_pes.xml
cp -f env_mach_pes.xml LockedFiles/env_mach_pes.xml.locked
cp -f env_build.xml LockedFiles/env_build.xml.locked

# Go back to original default layout
./cesm_setup -clean
cp env_build.xml.prs1 env_build.xml
cp env_mach_pes.xml.prs1 env_mach_pes.xml
./cesm_setup

