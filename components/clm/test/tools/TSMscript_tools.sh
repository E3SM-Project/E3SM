#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TSMscript_tools.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMscript_tools.$1.$2.$3.$4

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMscript_tools.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMscript_tools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMscript_tools.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMscript_tools.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSMscript_tools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

# Copy any sample files so can use them
cp $cfgdir/sample_* $rundir

optfile=${4%^*}
cfgfile=${4#*^}

if [[ "$2" == "PTCLM" ]]; then
  echo "TSMscript_tools.sh: calling TCBscripttools.sh to prepare executables for $2"
  ${CLM_SCRIPTDIR}/TCBscripttools.sh $1 $2 $cfgfile
  rc=$?
  if [ $rc -ne 0 ]; then
      echo "TSMscript_tools.sh: error from TCBscripttools.sh= $rc"
      echo "FAIL.job${JOBID}" > TestStatus
      exit 4
  fi 
  # Copy map files so we can use them
  subdir=1x1pt_US-UMB
  mkdir $rundir/$subdir
  cp $CSMDATA/lnd/clm2/PTCLMmydatafiles/$subdir/map_* $rundir/$subdir
elif [ "$optfile" != "$4" ]; then
  echo "TSMscript_tools.sh: calling TCBtools.sh to prepare $1 $2 executable"
  ${CLM_SCRIPTDIR}/TCBtools.sh $1 $2 $cfgfile
  rc=$?
  if [ $rc -ne 0 ]; then
      echo "TSMscript_tools.sh: error from TCBtools.sh= $rc"
      echo "FAIL.job${JOBID}" > TestStatus
      exit 4
  fi 
  tcbtools=${CLM_TESTDIR}/TCBtools.$1.$2.$cfgfile
else
  tcbtools="."
fi

scopts=`cat ${CLM_SCRIPTDIR}/nl_files/$optfile | sed -e "s|CSMDATA|$CSMDATA|g" | sed -e "s|EXEDIR|$tcbtools|" | sed -e "s|CFGDIR|$cfgdir|g"`

echo "TSMscript_tools.sh: running ${cfgdir}/$3 with $scopts; output in ${rundir}/test.log" 

if [ ! -f "${cfgdir}/$3" ]; then
    echo "TSMscript_tools.sh: error ${cfgdir}/$3 input script not found"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
   ${cfgdir}/$3 $scopts >> test.log 2>&1
   status="PASS"
   rc=$?
else
   echo "success" > test.log
   status="GEN"
   rc=0
fi

if [ $rc -eq 0 ] && grep -ci "success" test.log > /dev/null; then
    echo "TSMscript_tools.sh: smoke test passed" 
    echo "$status" > TestStatus
    # Copy files from subdirectories up...
    # (use hard links rather than symbolic links because 'ln -s' does funny
    # things when there are no matching files)
    ln */*.nc */*/*.nc .
else
    echo "TSMscript_tools.sh: error running $3, error= $rc" 
    echo "TSMscript_tools.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

exit 0
