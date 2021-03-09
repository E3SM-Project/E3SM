#!/bin/sh 
#

if [ -z "$BL_ROOT" ] && [ -z "$BL_TESTDIR" ]; then
    echo "TPF.sh: no environment variables set for baseline performance test - will skip" 
    exit 255
fi

if [ $# -ne 4 ]; then
    echo "TPF.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TPF.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TPF.sh: baseline performance test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}"
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TPF.sh: baseline performance test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TPF.sh: this baseline performance test failed under job ${prev_jobid} - moving those results to "
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
    echo "TPF.sh: error, unable to create work subdirectory" 
    exit 3
fi

cd ${rundir}

echo "TPF.sh: calling TSM.sh for smoke test"
${CAM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TPF.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ -n "${BL_ROOT}" ]; then
    if [ -z "$BL_TESTDIR" ]; then
	BL_TESTDIR=${CAM_TESTDIR}.bl
    fi
    echo "TPF.sh: generating baseline data from root $BL_ROOT - results in $BL_TESTDIR"

    echo "TPF.sh: calling ****baseline**** TSM.sh for smoke test"
    env CAM_TESTDIR=${BL_TESTDIR} \
	CAM_SCRIPTDIR=${BL_ROOT}/models/atm/cam/test/system \
	${BL_ROOT}/models/atm/cam/test/system/TSM.sh $1 $2 $3 $4
    rc=$?
    if [ $rc -ne 0 ]; then
	echo "TPF.sh: error from *baseline* TSM.sh= $rc" 
	echo "FAIL.job${JOBID}" > TestStatus
	exit 5
    fi
fi

if [ $4 = "build_only" ]; then
  exit 0
fi

#made the following extraction of timing info compatible with single or multiple timing files
model_time=`grep -i DRIVER_RUN_LOOP ${CAM_TESTDIR}/TSM.$1.$2.$3/ccsm_timing* | perl -e 'while (my $ll = <>) \
    { if ($ll =~ /DRIVER_RUN_LOOP[\s-]+[0-9]+[\s-]+([0-9\.]+)/) \
    { print "$1"; last }}' `
if [ "$model_time" = "" ]; then
    echo "TPF.sh: unable to determine model time from timing file(s) in ${CAM_TESTDIR}/TSM.$1.$2.$3/" 
    exit 6
fi

echo "TPF.sh: model_time= $model_time"

#made the following extraction of timing info compatible with single or multiple timing files
bl_time=`grep -i DRIVER_RUN_LOOP ${BL_TESTDIR}/TSM.$1.$2.$3/ccsm_timing* | perl -e 'while (my $ll = <>) \
    { if ($ll =~ /DRIVER_RUN_LOOP[\s-]+[0-9]+[\s-]+([0-9\.]+)/) \
    { print "$1"; last }}' `
if [ "$bl_time" = "" ]; then
    echo "TPF.sh: unable to determine model time from timing file(s) in ${BL_TESTDIR}/TSM.$1.$2.$3/" 
    exit 7
fi

echo "TPF.sh: bl_time= $bl_time"

perf_hit=`echo 5 k $model_time $bl_time / 100.00 * p | dc`
perf_hit=${perf_hit%.*}
perf_hit=`echo $perf_hit 100 - p | dc`
echo "TPF.sh: baseline performance hit= ${perf_hit}%   (cutoff is 5%)" 
if [ $perf_hit -gt 5 ]; then
    echo "TPF.sh: baseline performance test failed"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
else
    echo "TPF.sh: baseline performance test passed" 
    echo "PASS" > TestStatus
fi

exit 0
