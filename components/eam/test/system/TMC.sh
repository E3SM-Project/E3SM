#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TMC.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TMC.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TMC.sh: mass conservation test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}"
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TMC.sh: mass conservation test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TMC.sh: this mass conservation test failed under job ${prev_jobid} - moving those results to "
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
    echo "TMC.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TMC.sh: calling TSM.sh for smoke test"
${CAM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TMC.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ $4 = "build_only" ]; then
  exit 0
fi

result=`perl -e 'while (my $ll = <>) \
    { if ($ll =~ /before tphysbc DRY m=\s*[\d]+ name=TT_UN [^0-9]+([\S]+)/) \
    { print "$1 " }}' ${CAM_TESTDIR}/TSM.$1.$2.$3/test.log`

tstep=0
for next_val in $result; do
    if [ $tstep -eq 0 ]; then
	first_val=$next_val
    else
	if [ $next_val != $first_val ]; then
	    echo "TMC.sh: global dry TT_UN at step $tstep ( $next_val ) not equal to value at "
	    echo "                            step 0 ( $first_val )" 
	    echo "TMC.sh: see ${rundir}/test.log for details"
	    echo "FAIL.job${JOBID}" > TestStatus
	    exit 5
	fi
    fi
    tstep=`expr $tstep + 1`
done

if [ $tstep -eq 0 ]; then
    echo "TMC.sh: no mass conservation data available in output log"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
else
    echo "TMC.sh: mass conservation test passed" 
    echo "PASS" > TestStatus
fi

exit 0
