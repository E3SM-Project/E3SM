#!/bin/sh 
#

if [ -z "$BL_ROOT" ] && [ -z "$BL_TESTDIR" ]; then
    echo "TBL.sh: no environment variables set for baseline test - will skip" 
    exit 255
fi

if [ $# -ne 4 ]; then
    echo "TBL.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TBL.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBL.sh: baseline test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}"
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TBL.sh: baseline test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TBL.sh: this baseline test failed under job ${prev_jobid} - moving those results to "
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
    echo "TBL.sh: error, unable to create work subdirectory" 
    exit 3
fi

cd ${rundir}

echo "TBL.sh: calling TSM.sh for smoke test"
${CAM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBL.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ -n "${BL_ROOT}" ]; then
    if [ -z "$BL_TESTDIR" ]; then
	BL_TESTDIR=${CAM_TESTDIR}.bl
    fi
    echo "TBL.sh: generating baseline data from root $BL_ROOT - results in $BL_TESTDIR"

    echo "TBL.sh: calling ****baseline**** TSM.sh for smoke test"
    if [ "${CAM_BASEBACK}" = "YES" ]; then
        env CAM_TESTDIR=${BL_TESTDIR} \
  	    CAM_SCRIPTDIR=${BL_ROOT}/models/atm/cam/test/system \
	    ${BL_ROOT}/models/atm/cam/test/system/TSM.sh $1 $2 $3
    else
        env CAM_TESTDIR=${BL_TESTDIR} \
  	    CAM_SCRIPTDIR=${BL_ROOT}/models/atm/cam/test/system \
	    ${BL_ROOT}/models/atm/cam/test/system/TSM.sh $1 $2 $3 $4
    fi
    rc=$?
    if [ $rc -ne 0 ]; then
	echo "TBL.sh: error from *baseline* TSM.sh= $rc" 
	echo "FAIL.job${JOBID}" > TestStatus
	exit 5
    fi
fi

if [ $4 = "build_only" ]; then
  exit 0
fi

echo "TBL.sh: starting b4b comparisons " 
files_to_compare=`cd ${CAM_TESTDIR}/TSM.$1.$2.$3; ls *.cam*.h*.nc *.cam*.i*.nc`
if [ -z "${files_to_compare}" ]; then
    echo "TBL.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    ${CAM_SCRIPTDIR}/CAM_compare.sh \
	${BL_TESTDIR}/TSM.$1.$2.$3/${compare_file} \
	${CAM_TESTDIR}/TSM.$1.$2.$3/${compare_file}
    rc=$?
    mv cprnc.out cprnc.${compare_file}.out
    if [ $rc -eq 0 ]; then
        echo "TBL.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
	echo "TBL.sh: error from CAM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	all_comparisons_good="FALSE"
    fi
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TBL.sh: baseline test passed" 
    echo "PASS" > TestStatus
    if [ $CAM_RETAIN_FILES != "TRUE" ]; then
        echo "TBL.sh: removing some unneeded files to save disc space" 
        rm *.nc
        rm *.r*
    fi
else
    echo "TBL.sh: at least one file comparison did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

exit 0
