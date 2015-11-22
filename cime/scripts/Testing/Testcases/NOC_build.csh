#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here

./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_WAV  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_OCN  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_ICE  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_GLC  -val 2

set xmlquery_data=`./xmlquery -s JGFSEP CASE EXEROOT NTASKS_ATM NTASKS_LND NTASKS_ROF NTASKS_WAV NTASKS_OCN NTASKS_ICE NTASKS_GLC NTASKS_CPL`
set CASE=`echo $xmlquery_data | awk -F'JGFSEP' '{print $1}'`
set EXEROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $2}'`
set NTASKS_ATM=`echo $xmlquery_data | awk -F'JGFSEP' '{print $3}'`
set NTASKS_LND=`echo $xmlquery_data | awk -F'JGFSEP' '{print $4}'`
set NTASKS_ROF=`echo $xmlquery_data | awk -F'JGFSEP' '{print $5}'`
set NTASKS_WAV=`echo $xmlquery_data | awk -F'JGFSEP' '{print $6}'`
set NTASKS_OCN=`echo $xmlquery_data | awk -F'JGFSEP' '{print $7}'`
set NTASKS_ICE=`echo $xmlquery_data | awk -F'JGFSEP' '{print $8}'`
set NTASKS_GLC=`echo $xmlquery_data | awk -F'JGFSEP' '{print $9}'`
set NTASKS_CPL=`echo $xmlquery_data | awk -F'JGFSEP' '{print $10}'`

if ( $NTASKS_ATM == 1 ) then
  ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM  -val 2
endif
if ( $NTASKS_LND == 1 ) then
  ./xmlchange -file env_mach_pes.xml -id NTASKS_LND  -val 2
endif
if ( $NTASKS_ROF == 1 ) then
  ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF  -val 2
endif
if ( $NTASKS_WAV == 1 ) then
  ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV  -val 2
endif
if ( $NTASKS_ICE == 1 ) then
  ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE  -val 2
endif
if ( $NTASKS_GLC == 1 ) then
  ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC  -val 2
endif

./cesm_setup -clean
./cesm_setup

./$CASE.build
if ($status != 0) then
   exit -1    
endif 

