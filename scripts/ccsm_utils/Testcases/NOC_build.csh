#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

cd $CASEROOT

./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_WAV  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_OCN  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_ICE  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_GLC  -val 2
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

