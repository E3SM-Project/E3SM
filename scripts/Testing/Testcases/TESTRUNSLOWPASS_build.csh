#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

set EXEROOT  = `./xmlquery EXEROOT -value`
set RUNDIR   = `./xmlquery RUNDIR -value`
set TESTID   = `./xmlquery TEST_TESTID -value`

echo "#! /bin/bash" > $EXEROOT/${CIME_MODEL}.exe
echo "sleep 300" >> $EXEROOT/${CIME_MODEL}.exe
echo "echo Slow pass" >> $EXEROOT/${CIME_MODEL}.exe
echo "echo SUCCESSFUL TERMINATION > $RUNDIR/cpl.log.$TESTID" >> $EXEROOT/${CIME_MODEL}.exe

chmod +x $EXEROOT/${CIME_MODEL}.exe

echo "Build phase complete, just made simple script for ${CIME_MODEL}.exe" | tee $EXEROOT/cesm.bldlog.$TESTID

./xmlchange -noecho -file env_build.xml -id BUILD_COMPLETE -val TRUE

rm -rf LockedFiles/*

exit 0
