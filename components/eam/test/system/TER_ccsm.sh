#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TER_ccsm.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TER_ccsm.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TER_ccsm.sh: CESM exact restart test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TER_ccsm.sh: CESM exact restart test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TER_ccsm.sh: this CESM exact restart test failed under job ${prev_jobid} - moving those results to "
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
    echo "TER_ccsm.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

initial_length=${3%+*}
restart_string=${3#*+}
if [ ${initial_length} = $3 ] || [ ${restart_string} = $3 ]; then
    echo "TER_ccsm.sh: error processing input argument for run lengths"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

stop_flag=${restart_string##*[0-9]}
restart_length=${restart_string%%[sdm]}
if [ ${stop_flag} = ${restart_string} ] || [ ${restart_length} = ${restart_string} ]; then
    echo "TER_ccsm.sh: error processing input argument for run length= $3"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

case $stop_flag in
    s )  stop_option="nsteps";;

    d )  stop_option="ndays";;

    m )  stop_option="nmonths";;
esac

full_length=`expr $initial_length + $restart_length`

echo "TER_ccsm.sh: calling TSM_ccsm.sh for smoke test of full length ${full_length}${stop_flag}" 
${CAM_SCRIPTDIR}/TSM_ccsm.sh $1 $2 "${full_length}${stop_flag}" $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER_ccsm.sh: error from TSM_ccsm.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

if [ $4 = "build_only" ]; then
  exit 0
fi

echo "TER_ccsm.sh: calling TSM_ccsm.sh for smoke test of initial length ${initial_length}${stop_flag}" 
${CAM_SCRIPTDIR}/TSM_ccsm.sh $1 $2 "${initial_length}${stop_flag}" $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER_ccsm.sh: error from TSM_ccsm.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

echo "TER_ccsm.sh: restarting CESM; output in ${CAM_TESTDIR}/${test_name}/test.log" 

cd ${CAM_TESTDIR}/case.$1.$2
./xmlchange -file env_run.xml -id STOP_N -val $restart_length -silent
./xmlchange -file env_run.xml -id CONTINUE_RUN -val TRUE -silent
runscript=`ls *.run`
./$runscript > ${CAM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
cd ${rundir}
if [ $rc -eq 0 ] && grep -c "SUCCESSFUL TERMINATION" test.log > /dev/null; then
    echo "TER_ccsm.sh: restart of CESM completed successfully"
else
    echo "TER_ccsm.sh: error on restart run of CESM, error= $rc" 
    echo "TER_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

cp ${CAM_TESTDIR}/case.$1.$2/run/*cam.h0* .

echo "TER_ccsm.sh: starting b4b comparisons " 
files_to_compare=`ls *.cam*.h*.nc *.cam*.i*.nc`
if [ -z "${files_to_compare}" ]; then
    echo "TER_ccsm.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 9
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    ${CAM_SCRIPTDIR}/CAM_compare.sh \
	${compare_file} \
	${CAM_TESTDIR}/TSM_ccsm.$1.$2.${full_length}${stop_flag}/${compare_file}
    rc=$?
    mv cprnc.out cprnc.${compare_file}.out
    if [ $rc -eq 0 ]; then
        echo "TER_ccsm.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
	echo "TER_ccsm.sh: error from CAM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	all_comparisons_good="FALSE"
    fi
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TER_ccsm.sh: exact restart test passed" 
    echo "PASS" > TestStatus
    if [ $CAM_RETAIN_FILES != "TRUE" ]; then
        echo "TER_ccsm.sh: removing some unneeded files to save disc space" 
	rm *.nc
    fi
else
    echo "TER_ccsm.sh: at least one file comparison did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10
fi

exit 0
