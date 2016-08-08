#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

set EXEROOT  = `./xmlquery EXEROOT -value`
set RUNDIR   = `./xmlquery RUNDIR -value`
set TESTID   = `./xmlquery TEST_TESTID -value`

echo "#! /bin/bash" > $EXEROOT/cesm.exe
echo "sleep 300" >> $EXEROOT/cesm.exe
echo "echo Slow pass" >> $EXEROOT/cesm.exe
echo "echo SUCCESSFUL TERMINATION > $RUNDIR/cpl.log.$TESTID" >> $EXEROOT/cesm.exe

chmod +x $EXEROOT/cesm.exe

echo "Build phase complete, just made simple script for cesm.exe" | tee $EXEROOT/cesm.bldlog.$TESTID

./xmlchange -noecho -file env_build.xml -id BUILD_COMPLETE -val TRUE

rm -rf LockedFiles/*

exit 0
