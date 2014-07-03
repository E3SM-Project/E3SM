#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TBLscript_tools.sh: incorrect number of input arguments" 
    exit 1
fi

if [ -z "$BL_ROOT" ] && [ -z "$BL_TESTDIR" ]; then
    echo "TBLscript_tools.sh: no environment variables set for baseline test - will skip"
    exit 255
fi

test_name=TBLscript_tools.$1.$2.$3.$4

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBLscript_tools.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBLscript_tools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TBLscript_tools.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TBLscript_tools.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

rundir=${CLM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -r ${rundir}
fi
mkdir -p ${rundir} 
if [ $? -ne 0 ]; then
    echo "TBLscript_tools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TBLscript_tools.sh: calling TSMscript_tools.sh to run $1 $2 executable" 
${CLM_SCRIPTDIR}/TSMscript_tools.sh $1 $2 $3 $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBLscript_tools.sh: error from TSMtools.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ -n "${BL_ROOT}" ]; then
    if [ -z "$BL_TESTDIR" ]; then
        BL_TESTDIR=${CLM_TESTDIR}.bl
    fi
    echo "TBLscript_tools.sh: generating baseline data from root $BL_ROOT - results in $BL_TESTDIR"

    echo "TBLscript_tools.sh: calling ****baseline**** TSMtools.sh for smoke test"
    bl_dir=`/bin/ls -1d ${BL_ROOT}/models/lnd/clm/test/tools`
    env CLM_TESTDIR=${BL_TESTDIR} \
        CLM_SCRIPTDIR=$bl_dir \
        CLM_ROOT=$BL_ROOT \
        $bl_dir/TSMscript_tools.sh $1 $2 $3 $4
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "TBLscript_tools.sh: error from *baseline* TSMscript_tools.sh= $rc"
        echo "FAIL.job${JOBID}" > TestStatus
        exit 5
    fi
fi

echo "TBLscript_tools.sh: starting b4b comparisons "
files_to_compare=`cd ${CLM_TESTDIR}/TSMscript_tools.$1.$2.$3.$4; ls *.nc`
if [ -z "${files_to_compare}" ] && [ "$debug" != "YES" ]; then
    echo "TBLscript_tools.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    env CPRNC_EXE=${CLM_SCRIPTDIR}/../../tools/shared/ncl_scripts/cprnc.pl \
        ${CLM_SCRIPTDIR}/CLM_compare.sh \
        ${BL_TESTDIR}/TSMscript_tools.$1.$2.$3.$4/${compare_file} \
        ${CLM_TESTDIR}/TSMscript_tools.$1.$2.$3.$4/${compare_file}
    rc=$?
    mv cprnc.out cprnc.${compare_file}.out
    if [ $rc -eq 0 ]; then
        echo "TBLscript_tools.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
        echo "TBLscript_tools.sh: error from CLM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details
"
        all_comparisons_good="FALSE"
    fi
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TBLscript_tools.sh: baseline test passed"
    echo "PASS" > TestStatus
    if [ $CLM_RETAIN_FILES != "TRUE" ]; then
        echo "TBLscript_tools.sh: removing some unneeded files to save disc space"
        rm *.nc
        rm *.r*
    fi
else
    echo "TBLscript_tools.sh: at least one file comparison did not pass"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi



exit 0
