#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

set xmlquery_data=`./xmlquery -s JGFSEP EXEROOT TEST_TESTID`
set EXEROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $1}'`
set TESTID=`echo $xmlquery_data | awk -F'JGFSEP' '{print $2}'`

echo "ERROR: Intentional fail for testing infrastructure" | tee $EXEROOT/cesm.bldlog.$TESTID

exit -1
