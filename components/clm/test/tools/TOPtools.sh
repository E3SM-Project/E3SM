#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TOPtools.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TOPtools.$1.$2.$3.$4

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TOPtools.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TOPtools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TOPtools.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TOPtools.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
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
    echo "TOPtools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

if [ ${CLM_THREADS} -lt 2 ]; then
   echo "TOPtools.sh: error not enough threads are being used to do the comparision"
   echo "FAIL.job${JOBID}" > TestStatus
   exit 5
fi
if [ "$3" != "tools__o" ] && [ "$3" != "tools__do" ]; then
   echo "TOPtools.sh: error build needs to be done Open-MP"
   echo "FAIL.job${JOBID}" > TestStatus
   exit 5
fi

echo "TOPtools.sh: calling TSMtools.sh to run $1 $2 executable"
${CLM_SCRIPTDIR}/TSMtools.sh $1 $2 $3 $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TOPtools.sh: error from TSMtools.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi
mkdir $rundir/$CLM_THREADS
cp ${CLM_TESTDIR}/TSMtools.$1.$2.$3.$4/*.nc $rundir/$CLM_THREADS

# Get a list of different threads to run for, powers of 2 from 1 up to the thread count
threads=1
list="1 "
until [ "$threads" -ge "$CLM_THREADS" ]; do 
  threads=`perl -e "$CLM_THREADS<$threads*2 ? print $CLM_THREADS : print $threads*2"`
  if [ "$threads" -lt "$CLM_THREADS" ]; then list="$list $threads "; fi
done

all_comparisons_good="TRUE"
for threads in $list
do 
  echo "TOPtools.sh: calling TSMtools.sh to run $1 executable for $threads threads" 
  env CLM_THREADS=$threads CLM_RERUN=yes ${CLM_SCRIPTDIR}/TSMtools.sh $1 $2 $3 $4
  rc=$?
  if [ $rc -ne 0 ]; then
      echo "TOPtools.sh: error from TSMtools.sh= $rc" 
      echo "FAIL.job${JOBID}" > TestStatus
      exit 6
  fi
  mkdir $rundir/$threads
  cp ${CLM_TESTDIR}/TSMtools.$1.$2.$3.$4/*.nc $rundir/$threads
  files_to_compare=`cd $rundir/$threads; ls *.nc`
  for compare_file in ${files_to_compare}; do

      env CPRNC_EXE=${CLM_SCRIPTDIR}/../../tools/shared/ncl_scripts/cprnc.pl \
          ${CLM_SCRIPTDIR}/CLM_compare.sh \
          $rundir/$CLM_THREADS/${compare_file} \
          $rundir/$threads/${compare_file}
      rc=$?
      cprout="cprnc.${compare_file}.threads${threads}.out"
      mv cprnc.out $cprout
      if [ $rc -eq 0 ]; then
          echo "TOPtools.sh: comparison successful; output in $cprout"
      else
          echo "TOPtools.sh: error from CLM_compare.sh= $rc; see $cprout for details"
          all_comparisons_good="FALSE"
      fi
  done
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TOPtools.sh: OpenMP comparison test passed"
    echo "PASS" > TestStatus
    if [ $CLM_RETAIN_FILES != "TRUE" ]; then
        echo "TOPtools.sh: removing some unneeded files to save disc space"
        rm */*.nc
    fi
else
    echo "TOPtools.sh: at least one file comparison did not pass"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

exit 0
