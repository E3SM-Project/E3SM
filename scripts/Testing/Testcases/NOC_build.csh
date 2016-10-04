#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE        = `./xmlquery CASE		-value`
set EXEROOT     = `./xmlquery EXEROOT		-value`

./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_WAV  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_OCN  -val 1
./xmlchange -file env_mach_pes.xml -id NINST_ICE  -val 2
./xmlchange -file env_mach_pes.xml -id NINST_GLC  -val 2

set NTASKS_ATM  = `./xmlquery NTASKS_ATM	-value`
set NTASKS_LND  = `./xmlquery NTASKS_LND	-value`
set NTASKS_ROF  = `./xmlquery NTASKS_ROF	-value`
set NTASKS_WAV  = `./xmlquery NTASKS_WAV	-value`
set NTASKS_OCN  = `./xmlquery NTASKS_OCN	-value`
set NTASKS_ICE  = `./xmlquery NTASKS_ICE	-value`
set NTASKS_GLC  = `./xmlquery NTASKS_GLC	-value`
set NTASKS_CPL  = `./xmlquery NTASKS_CPL	-value`

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

