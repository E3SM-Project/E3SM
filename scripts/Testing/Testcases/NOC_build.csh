#!/bin/csh -f
setenv CIMEROOT `./xmlquery CIMEROOT    --value`
if( -e env_mach_pes.xml.orig) then
    cp env_mach_pes.xml.orig env_mach_pes.xml
    ./case.setup -clean
else
    cp env_mach_pes.xml env_mach_pes.xml.orig
endif
./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE        = `./xmlquery CASE		--value`
set EXEROOT     = `./xmlquery EXEROOT		--value`

./xmlchange --file env_mach_pes.xml --id NINST_OCN  --val 1
set maxtasks = 0

foreach comp (ATM LND ROF WAV ICE GLC)
  ./xmlchange NINST_$comp=2

  set tasks = `./xmlquery NTASKS_$comp --value`
  set root = `./xmlquery ROOTPE_$comp --value`
  set even = `expr $tasks % 2`
  if( $even == 1) then
    set tasks = `expr $tasks \* 2`
    ./xmlchange NTASKS_$comp=$tasks
  endif
  if( $root > 0) then
    set rooteven = `expr $root % 2`
    if( $rooteven == 1 ) then
      set root = `expr $root \* 2`
      ./xmlchange ROOTPE_$comp=$root
    endif
    set mtasks = `expr $root + $tasks`
    if ( $mtasks > $maxtasks ) then
      set maxtasks = $mtasks
    endif
  endif
end
set ocnroot = `./xmlquery ROOTPE_OCN --value`
if( ( $ocnroot > 0 ) && ( $ocnroot < $maxtasks ) ) then
  ./xmlchange ROOTPE_OCN=$maxtasks
endif

./case.setup -clean
./case.setup

./case.build --testmode
if ($status != 0) then
   exit -1    
endif 

