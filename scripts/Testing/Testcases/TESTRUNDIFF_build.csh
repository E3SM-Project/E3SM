#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

set EXEROOT  = `./xmlquery EXEROOT --value`
set RUNDIR   = `./xmlquery RUNDIR --value`
set TESTID   = `./xmlquery TEST_TESTID --value`
set CIMEROOT = `./xmlquery CIMEROOT --value`
set CASE     = `./xmlquery CASE --value`

echo "#! /bin/bash" > $EXEROOT/${CIME_MODEL}.exe
echo "echo Insta pass" >> $EXEROOT/${CIME_MODEL}.exe
echo -n "echo SUCCESSFUL TERMINATION > $RUNDIR/cpl.log." >> $EXEROOT/${CIME_MODEL}.exe
echo '$LID' >> $EXEROOT/${CIME_MODEL}.exe
echo "cp $CIMEROOT/utils/python/tests/cpl.hi1.nc.test $RUNDIR/$CASE.cpl.hi.0.nc.base" >> $EXEROOT/${CIME_MODEL}.exe

chmod +x $EXEROOT/${CIME_MODEL}.exe

echo "Build phase complete, just made simple script for ${CIME_MODEL}.exe" | tee $EXEROOT/cesm.bldlog.$TESTID

./xmlchange --noecho --file env_build.xml --id BUILD_COMPLETE --val TRUE

rm -rf LockedFiles/*

exit 0
