#!/bin/csh -f
setenv CIMEROOT `./xmlquery CIMEROOT    --value`

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE     = `./xmlquery CASE    --value`

./case.build $*
if ($status != 0) then
   exit -1
endif

