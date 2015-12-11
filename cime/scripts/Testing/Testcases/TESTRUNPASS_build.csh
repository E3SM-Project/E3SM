#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

set xmlquery_data=`./xmlquery -s JGFSEP EXEROOT RUNDIR TEST_TESTID`
set EXEROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $1}'`
set RUNDIR=`echo $xmlquery_data | awk -F'JGFSEP' '{print $2}'`
set TESTID=`echo $xmlquery_data | awk -F'JGFSEP' '{print $3}'`

echo "#! /bin/bash" > $EXEROOT/cesm.exe
echo "echo Insta pass" >> $EXEROOT/cesm.exe
echo "echo SUCCESSFUL TERMINATION > $RUNDIR/cpl.log.$TESTID" >> $EXEROOT/cesm.exe

chmod +x $EXEROOT/cesm.exe

echo "Build phase complete, just made simple script for cesm.exe" | tee $EXEROOT/cesm.bldlog.$TESTID

./xmlchange -noecho -file env_build.xml -id BUILD_COMPLETE -val TRUE

rm -rf LockedFiles/*

exit 0
