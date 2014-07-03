#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TSMtools.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMtools.$1.$2.$3.$4

if [ -z "$CLM_RERUN" ]; then
  CLM_RERUN="no"
fi

if [ "$CLM_RERUN" != "yes" ] && [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMtools.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMtools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMtools.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMtools.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=`ls -1d ${CLM_ROOT}/models/lnd/clm/tools/$1/$2`
rundir=${CLM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -r ${rundir}
fi
mkdir -p ${rundir} 
if [ $? -ne 0 ]; then
    echo "TSMtools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "Copy any text files over"
cp $cfgdir/*.txt $rundir

echo "TSMtools.sh: calling TCBtools.sh to prepare $1 $2 executable" 
${CLM_SCRIPTDIR}/TCBtools.sh $1 $2 $3
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMtools.sh: error from TCBtools.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

echo "TSMtools.sh: running $1; output in ${rundir}/test.log" 

if [ "$3" = "tools__o" ] || [ "$3" = "tools__do" ]; then
   toolrun="env OMP_NUM_THREADS=${CLM_THREADS} ${CLM_TESTDIR}/TCBtools.$1.$2.$3/$2"
else
   toolrun="${CLM_TESTDIR}/TCBtools.$1.$2.$3/$2"
fi

runfile="${cfgdir}/$2.$4"

if [ ! -f "${runfile}" ]; then
    runfile="${CLM_SCRIPTDIR}/nl_files/$2.$4"
    if [ ! -f "${runfile}" ]; then
        echo "TSMtools.sh: error ${runfile} input run file not found"
        echo "FAIL.job${JOBID}" > TestStatus
        exit 5
    fi
fi

echo "Run file type = ${4#*.}"
if [ ${4#*.} == "runoptions" ]; then
  echo "$toolrun "`cat ${runfile}`
  cp $cfgdir/*.nc .
  if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
     $toolrun  `cat ${runfile}` >> test.log 2>&1
     rc=$?
     status="PASS"
  else
     echo "Successfully created file" > test.log
     status="GEN"
     rc=0
  fi
else
  echo "$toolrun < ${runfile}"
  if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
     $toolrun < ${runfile} >> test.log 2>&1
     rc=$?
     status="PASS"
  else
     echo "Successfully created file" > test.log
     status="GEN"
     rc=0
  fi
fi

if [ $rc -eq 0 ] && grep -ci "Successfully created " test.log > /dev/null; then
    echo "TSMtools.sh: smoke test passed" 
    echo "$status" > TestStatus
else
    echo "TSMtools.sh: error running $1 $2, error= $rc" 
    echo "TSMtools.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

exit 0
