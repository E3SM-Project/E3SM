#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TCBtools.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCBtools.$1.$2.$3

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCBtools.sh: build test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCBtools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCBtools.sh: build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCBtools.sh: this build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=`ls -1d ${CLM_ROOT}/models/lnd/clm/tools/$1/$2`
blddir=${CLM_TESTDIR}/${test_name}/src
if [ -d ${blddir} ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir}
if [ $? -ne 0 ]; then
    echo "TCBtools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

echo "TCBtools.sh: building $1 $2 executable; output in ${blddir}/test.log" 
#
# Copy build files over
#
cp $cfgdir/src/Makefile      .
cp $cfgdir/src/Srcfiles      .
cp $cfgdir/src/Mkdepends     .
cp $cfgdir/src/Makefile.common .
#
# Add cfgdir path to beginning of each path in Filepath
#
touch Filepath
while read filepath_arg; do
    echo "${cfgdir}/src/${filepath_arg}" >> Filepath
done < ${cfgdir}/src/Filepath

#
# Figure out configuration
#
if [ ! -f ${CLM_SCRIPTDIR}/config_files/$3 ]; then
    echo "TCB.sh: configure options file ${CLM_SCRIPTDIR}/config_files/$3 not found"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

##construct string of args to configure
config_string="$TOOLS_MAKE_STRING TOOLROOT=$cfgdir "
while read config_arg; do
    config_string="${config_string}${config_arg} "
done < ${CLM_SCRIPTDIR}/config_files/$3

attempt=1
still_compiling="TRUE"
if [ "$TOOLSLIBS" != "" ]; then
   export SLIBS=$TOOLSLIBS
fi
while [ $still_compiling = "TRUE" ]; do

   if [ "$2" = "gen_domain" ]; then
      HOSTNAME=`uname -n | cut -c 1-2`
      if [ "$HOSTNAME" = "be" ]; then
         echo "TCBtools.sh: run configure for gen_domain on bluefire"
         env CCSMROOT=${CLM_ROOT} ${CLM_ROOT}/scripts/ccsm_utils/Machines/configure -mach bluefire >> test.log 2>&1
         rc=$?
      fi
   fi

    echo "TCBtools.sh: call to make:" 
    echo "        ${MAKE_CMD} ${config_string} "
    if [ "$debug" != "YES" ]; then
       ${MAKE_CMD} ${config_string} >> test.log 2>&1
       status="PASS"
       rc=$(( $rc + $? ))
    else
       status="GEN"
       rc=0
    fi
    if [ $rc -eq 0 ]; then
	echo "TCBtools.sh: make was successful" 
	echo "TCBtools.sh: configure and build test passed"
	echo "$status" > TestStatus
	if [ $CLM_RETAIN_FILES != "TRUE" ]; then
	    echo "TCBtools.sh: removing some unneeded files to save disc space" 
	    rm *.o
	    rm *.mod
	fi
	still_compiling="FALSE"
    elif [ $attempt -lt 10 ] && \
        grep -c "LICENSE MANAGER PROBLEM" test.log > /dev/null; then
        attempt=`expr $attempt + 1`
        echo "TCBtools.sh: encountered License Manager Problem; launching attempt #$attempt"
    else
	echo "TCBtools.sh: clm build failed, error from make= $rc" 
	echo "TCBtools.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 6
    fi
done
if [ "$TOOLSLIBS" != "" ]; then
   export -n SLIBS
fi

exit 0
