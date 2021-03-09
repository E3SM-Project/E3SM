#!/bin/sh
#

if [ $# -ne 6 ]; then
    echo "TSC.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSC.$1.$2.$3.$4.$5

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSC.sh: scam b4b test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}"
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSC.sh: scam b4b test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSC.sh: this scam b4b test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSC.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSC.sh: calling TSM.sh for cam to generate iop datafiles"
${CAM_SCRIPTDIR}/TSM.sh $1 $2 $5 $6
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSC.sh: error from TSM.sh for cam run= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

#temporarily stage these files one level up in work directory tree
cp ${CAM_TESTDIR}/TSM.$1.$2.$5/camrun.*.i.*.nc ../.
cp ${CAM_TESTDIR}/TSM.$1.$2.$5/camrun.*.h1.*.nc ../.

echo "TSC.sh: calling TSM.sh using iop files as input to single-column model"
${CAM_SCRIPTDIR}/TSM.sh $3 $4 $5 $6
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSC.sh: error from TSM.sh for scam run= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi


if [ $6 = "build_only" ]; then
  exit 0
fi

#remove temporarily staged files
if [ $CAM_RETAIN_FILES != "TRUE" ]; then
    rm ../camrun.*
fi

# Now test the output
echo "TSC.sh: Comparing answers to ensure SCAM gives bit-for-bit answers as CAM ... "
myvar=`ncdump -ff -p 9,17 -v QDIFF,TDIFF ${CAM_TESTDIR}/TSM.$3.$4.$5/camrun.*.h1.*00000.nc | egrep //\.\*DIFF | sed s/^\ \*// | sed s/\[,\;\].\*\$// | uniq`
if [ "$myvar" == "0" ]; then
    myvar=`ncdump -ff -p 9,17 -v QDIFF,TDIFF ${CAM_TESTDIR}/TSM.$3.$4.$5/camrun.*.h1.*08400.nc | egrep //\.\*DIFF | sed s/^\ \*// | sed s/\[,\;\].\*\$// | uniq`
    if [ "$myvar" == "0" ]; then
      echo "TSC.sh:  scam b4b test passed" 
      echo "PASS" > TestStatus
    else
      echo "TSC.sh: scam b4b test did not pass" 
      echo "FAIL.job${JOBID}" > TestStatus
      exit 6
    fi
else
    echo "TSC.sh: scam b4b test did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

exit 0
