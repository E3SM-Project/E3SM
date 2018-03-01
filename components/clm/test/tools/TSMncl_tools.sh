#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "TSMncl_tools.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMncl_tools.$1.$2

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMncl_tools.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMncl_tools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMncl_tools.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMncl_tools.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=`ls -1d ${CLM_ROOT}/models/lnd/clm/tools/$1/ncl_scripts`
rundir=${CLM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -r ${rundir}
fi
mkdir -p ${rundir} 
if [ $? -ne 0 ]; then
    echo "TSMncl_tools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSMncl_tools.sh: running $1 $2; output in ${rundir}/test.log" 

if [ ! -f "${cfgdir}/$2.ncl" ]; then
    echo "TSMncl_tools.sh: error ${cfgdir}/$2.ncl input script not found"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
   ncl ${cfgdir}/$1.ncl >> test.log 2>&1
   status="PASS"
   rc=$?
else
   echo "success" > test.log
   status="GEN"
   rc=0
fi

if [ $rc -eq 0 ] && grep -ci "success" test.log > /dev/null; then
    echo "TSMncl_tools.sh: smoke test passed" 
    echo "$status" > TestStatus
else
    echo "TSMncl_tools.sh: error running $1 $2, error= $rc" 
    echo "TSMncl_tools.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

exit 0
