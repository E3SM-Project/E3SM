#!/bin/csh -f

# Need the following environment variables to be available
# $CCSMROOT $CASEROOT $CASE $TESTCASE
# $GENERATE_BASELINE, $COMPARE_BASELINE, $STOP_OPTION, $STOP_N

set currdir = `pwd`

cd ${CASEROOT} || exit -1
source ./Tools/ccsm_getenv || exit -2

#-------------------------------------------------------------
# Modify env_run.xml in test directory 
#-------------------------------------------------------------
if ( $?GENERATE_BASELINE || $?COMPARE_BASELINE ) then
  ./xmlchange -file env_run.xml -id HIST_OPTION -val '$STOP_OPTION'
  ./xmlchange -file env_run.xml -id HIST_N      -val '$STOP_N'
endif

#-------------------------------------------------------------
# create run script (that testing script will call)
#-------------------------------------------------------------
./cesm_setup
if ($status != 0) then
  echo "ERROR: testcase_setup.csh invocation of cesm_setup failed"
  exit -1
endif

#-------------------------------------------------------------
# create testing script
#-------------------------------------------------------------
cd $CASEROOT || exit -2

# create batch header for testing script
env PHASE=set_batch env TESTMODE=test ${CCSM_MACHDIR}/mkbatch.${MACH} || exit -3  

# set environment variables for testing...
source ${CCSMROOT}/scripts/ccsm_utils/Tools/testcase_env.csh

# fill in the rest of the test script...
cat >> ${CASE}.test << EOF

cd $CASEROOT
./Tools/ccsm_check_lockedfiles || exit -4
source ./Tools/ccsm_getenv || exit -5

EOF

cat ${CCSMROOT}/scripts/ccsm_utils/Tools/testcase_begin         >> $CASE.test || exit -5
cat ${CCSMROOT}/scripts/ccsm_utils/Testcases/${TESTCASE}_script >> $CASE.test || exit -6
cat ${CCSMROOT}/scripts/ccsm_utils/Tools/testcase_end           >> $CASE.test || exit -7
if (-e ${CCSMROOT}/scripts/ccsm_utils/Testcases/${TESTCASE}_build.csh) then
  cp -f ${CCSMROOT}/scripts/ccsm_utils/Testcases/${TESTCASE}_build.csh ./${CASE}.test_build
else
  cp -f ${CCSMROOT}/scripts/ccsm_utils/Testcases/tests_build.csh ./${CASE}.test_build
endif

# copy over compare script
#-------------------------------------------------------------
cp ${CCSMROOT}/scripts/ccsm_utils/Tools/check_exactrestart.pl ./Tools/ || exit -8
cp ${CCSMROOT}/scripts/ccsm_utils/Tools/check_memory.pl       ./Tools/ || exit -8
cp ${CCSMROOT}/scripts/ccsm_utils/Tools/compare_throughput.pl ./Tools/ || exit -8
cp ${CCSMROOT}/scripts/ccsm_utils/Tools/hist_compare.csh      ./Tools/ || exit -9

chmod 755 $CASE* *pl

cd $currdir

