#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

set EXEROOT  = `./xmlquery EXEROOT -value`
set TESTID   = `./xmlquery TEST_TESTID -value`

echo "ERROR: Intentional fail for testing infrastructure" | tee $EXEROOT/cesm.bldlog.$TESTID

exit -1
