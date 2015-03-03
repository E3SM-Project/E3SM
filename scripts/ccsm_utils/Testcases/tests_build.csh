#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE     = `./xmlquery CASE    -value`

./$CASE.build
if ($status != 0) then
   echo "Error: build failed" >! ./TestStatus
   echo "CFAIL $CASE" > ./TestStatus
   exit -1    
endif 

