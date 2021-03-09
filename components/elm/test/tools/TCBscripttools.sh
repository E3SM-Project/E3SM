#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TCBscripttools.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCBscripttools.$1.$2.$3

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCBscripttools.sh: build test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCBscripttools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCBscripttools.sh: build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCBscripttools.sh: this build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=`ls -1d ${CLM_ROOT}/models/lnd/clm/tools/$1/$2`
blddir=${CLM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir}
if [ $? -ne 0 ]; then
    echo "TCBscripttools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

echo "TCBscripttools.sh: building $2 executables; output in ${blddir}/test.log" 
#
# Build script to exercise
#
if [ ! -x ${cfgdir}/$3 ]; then
    echo "TCB.sh: build run script file ${cfgdir}/$3 not found"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

echo "TCBscripttools.sh: run the build scriptmake:" 
echo "        ${cfgdir}/$3"

if [ "$debug" != "YES" ]; then
    export CESM_ROOT=${CLM_ROOT}
    ${cfgdir}/$3  >> test.log 2>&1
    rc=$(( $rc + $? ))
    status="PASS"
else
    status="GEN"
    rc=0
fi
if [ $rc -eq 0 ]; then
    echo "TCBscripttools.sh: build script was successful" 
    echo "TCBscripttools.sh: build script test passed"
    echo "$status" > TestStatus
else
    echo "TCBscripttools.sh: clm build script failed, error from build script= $rc" 
    echo "TCBscripttools.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

exit 0
