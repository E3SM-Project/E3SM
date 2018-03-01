#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TSM_ccsm.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSM_ccsm.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSM_ccsm.sh: CESM smoke test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSM_ccsm.sh: CESM smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSM_ccsm.sh: this CESM smoke test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSM_ccsm.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSM_ccsm.sh: calling TCB_ccsm.sh to prepare cam executable" 
${CAM_SCRIPTDIR}/TCB_ccsm.sh $1 $2 $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSM_ccsm.sh: error from TCB_ccsm.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ $4 = "build_only" ]; then
  exit 0
fi

#if [ -f ${CAM_TESTDIR}/case.$1.$2/run/*cam.h0* ]; then
    rm -f ${CAM_TESTDIR}/case.$1.$2/run/*cam.h0*
#fi

stop_flag=${3##*[0-9]}
run_length=${3%%[sdm]}
if [ ${stop_flag} = $3 ] || [ ${run_length} = $3 ]; then
    echo "TSM_ccsm.sh: error processing input argument for run length= $3"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

case $stop_flag in
    s )  stop_option="nsteps";;

    d )  stop_option="ndays";;

    m )  stop_option="nmonths";;
esac

echo "TSM_ccsm.sh: running CESM; output in ${CAM_TESTDIR}/${test_name}/test.log" 

cd ${CAM_TESTDIR}/case.$1.$2
./xmlchange -file env_run.xml -id STOP_N -val $run_length -silent
./xmlchange -file env_run.xml -id STOP_OPTION -val $stop_option -silent
./xmlchange -file env_run.xml -id DOUT_S -val FALSE -silent
runscript=`ls *.run` 
./$runscript > ${CAM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
cd ${rundir}
if [ $rc -eq 0 ] && grep -c "SUCCESSFUL TERMINATION" test.log > /dev/null; then
    echo "TSM_ccsm.sh: CESM smoke test passed" 
    echo "PASS" > TestStatus
    if [ $CAM_RETAIN_FILES != "TRUE" ]; then
        echo "TSM_ccsm.sh: removing some unneeded files to save disc space" 
        #think of any?
    fi
else
    echo "TSM_ccsm.sh: error running CESM, error= $rc" 
    echo "TSM_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

cp ${CAM_TESTDIR}/case.$1.$2/run/*cam.h0* .

exit 0
