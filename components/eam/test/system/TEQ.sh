#!/bin/sh 
#

if [ $# -ne 6 ]; then
    echo "TEQ.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TEQ.$1.$2.$3.$4.$5

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TEQ.sh: equivalent run test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}"
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TEQ.sh: equivalent run test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TEQ.sh: this equivalent run test failed under job ${prev_jobid} - moving those results to "
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
    echo "TEQ.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TEQ.sh: calling TSM.sh for first smoke test"
${CAM_SCRIPTDIR}/TSM.sh $1 $2 $5 $6
rc=$?
if [ $rc -ne 0 ]; then
    echo "TEQ.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

echo "TEQ.sh: calling TSM.sh for second smoke test"
${CAM_SCRIPTDIR}/TSM.sh $3 $4 $5 $6
rc=$?
if [ $rc -ne 0 ]; then
    echo "TEQ.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

if [ $6 = "build_only" ]; then
  exit 0
fi

echo "TEQ.sh: starting b4b comparisons " 
files_to_compare=`cd ${CAM_TESTDIR}/TSM.$1.$2.$5; ls *.cam*.h*.nc *.cam*.i*.nc`
if [ -z "${files_to_compare}" ]; then
    echo "TEQ.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    ${CAM_SCRIPTDIR}/CAM_compare.sh \
	${CAM_TESTDIR}/TSM.$1.$2.$5/${compare_file} \
	${CAM_TESTDIR}/TSM.$3.$4.$5/${compare_file}
    rc=$?
    mv cprnc.out cprnc.${compare_file}.out
    if [ $rc -eq 0 ]; then
        echo "TEQ.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
	echo "TEQ.sh: error from CAM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	all_comparisons_good="FALSE"
    fi
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TEQ.sh: equivalent run test passed" 
    echo "PASS" > TestStatus
else
    echo "TEQ.sh: at least one file comparison did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

exit 0
