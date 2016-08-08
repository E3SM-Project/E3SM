#!/bin/csh -f

./Tools/check_lockedfiles || exit -1

setenv COMPILER `./xmlquery COMPILER -value`
source env_mach_specific || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set EXEROOT  = `./xmlquery EXEROOT -value`
set MACH     = `./xmlquery MACH -value`
set PROCS    = `./xmlquery TOTALPES -value`
set ACMEROOT = `./xmlquery CCSMROOT -value`
set BASELINE = `./xmlquery CCSM_BASELINE -value`
set BASEGEN  = `./xmlquery BASEGEN_CASE -value`
set BASECMP  = `./xmlquery BASECMP_CASE -value`
set TESTID   = `./xmlquery TEST_TESTID -value`
set GEN      = `./xmlquery GENERATE_BASELINE -value`

cd $EXEROOT

if ($GEN == "TRUE") then
    cmake -C $ACMEROOT/components/homme/cmake/machineFiles/${MACH}.cmake -DUSE_NUM_PROCS=$PROCS $ACMEROOT/components/homme -DHOMME_BASELINE_DIR=$BASELINE/$BASEGEN >& $EXEROOT/homme.bldlog || exit -1
else
    cmake -C $ACMEROOT/components/homme/cmake/machineFiles/${MACH}.cmake -DUSE_NUM_PROCS=$PROCS $ACMEROOT/components/homme -DHOMME_BASELINE_DIR=$BASELINE/$BASECMP >& $EXEROOT/homme.bldlog || exit -1
endif

make -j8 >>& $EXEROOT/homme.bldlog || exit -1

cd -

./xmlchange -noecho -file env_build.xml -id BUILD_COMPLETE -val TRUE

rm -rf LockedFiles/*

exit 0
