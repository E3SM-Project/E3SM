#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TDD.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TDD.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TDD.sh: Divergence damper test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}"
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TDD.sh: Divergence damper test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TDD.sh: this Divergence damper test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CAM_TESTDIR}/${test_name} ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

rundir=${CAM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -rf ${rundir}
fi
mkdir -p ${rundir} 
if [ $? -ne 0 ]; then
    echo "TDD.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TDD.sh: calling TSM.sh for smoke test"
${CAM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TDD.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ $4 = "build_only" ]; then
  exit 0
fi

if grep -c "Divergence damper for spectral dycore invoked"  ${CAM_TESTDIR}/TSM.$1.$2.$3/test.log > /dev/null; then
    echo "TDD.sh: Divergence Damper test passed" 
    echo "PASS" > TestStatus
else
    echo "TDD.sh: Divergence Damper test failed"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

exit 0
