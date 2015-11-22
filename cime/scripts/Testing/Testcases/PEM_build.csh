#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE        = `./xmlquery CASE`
set EXEROOT     = `./xmlquery EXEROOT`

# Reset all previous env_mach_pes settings

if ( -e env_mach_pes.xml.1 )  then
  cp -f env_mach_pes.xml.1 env_mach_pes.xml
endif

./cesm_setup -clean
./cesm_setup

./$CASE.build
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.1
cp -f env_build.xml      env_build.xml.1
cp -f env_mach_pes.xml   env_mach_pes.xml.1

set xmlquery_data=`./xmlquery -s JGFSEP NTASKS_ATM NTASKS_LND NTASKS_ROF NTASKS_WAV NTASKS_OCN NTASKS_ICE NTASKS_GLC NTASKS_CPL`
set NTASKS_ATM=`echo $xmlquery_data | awk -F'JGFSEP' '{print $1}'`
set NTASKS_LND=`echo $xmlquery_data | awk -F'JGFSEP' '{print $2}'`
set NTASKS_ROF=`echo $xmlquery_data | awk -F'JGFSEP' '{print $3}'`
set NTASKS_WAV=`echo $xmlquery_data | awk -F'JGFSEP' '{print $4}'`
set NTASKS_OCN=`echo $xmlquery_data | awk -F'JGFSEP' '{print $5}'`
set NTASKS_ICE=`echo $xmlquery_data | awk -F'JGFSEP' '{print $6}'`
set NTASKS_GLC=`echo $xmlquery_data | awk -F'JGFSEP' '{print $7}'`
set NTASKS_CPL=`echo $xmlquery_data | awk -F'JGFSEP' '{print $8}'`

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

./cesm_setup -clean -testmode
./cesm_setup

./$CASE.build
if ($status != 0) then
   exit -1    
endif 

mv -f $EXEROOT/cesm.exe $EXEROOT/cesm.exe.2
cp -f env_build.xml      env_build.xml.2
cp -f env_mach_pes.xml   env_mach_pes.xml.2



