#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TCB_ccsm.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCB_ccsm.$1.$2

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCB_ccsm.sh: CESM configure and build test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCB_ccsm.sh: CESM configure and build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCB_ccsm.sh: this CESM configure and build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CAM_TESTDIR}/${test_name} ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
        if [ $3 = "run_only" ]; then
            echo "TCB_ccsm.sh: CAM_RBOPTIONS set to run_only, this test needs to be built.  Try build_only or run_and_build"
            exit 2
        fi
    fi
fi

blddir=${CAM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -rf ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCB_ccsm.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

echo "TCB_ccsm.sh: building ccsm executable; output in ${CAM_TESTDIR}/${test_name}/test.log" 
if [ -d ${CAM_TESTDIR}/case.$1.$2 ]; then
    rm -rf ${CAM_TESTDIR}/case.$1.$2
fi

# determine if chemistry preprocessor needs to be invoked
compset=${2%+*}
usrmech=${2#*+}

${CAM_ROOT}/scripts/create_newcase -case ${CAM_TESTDIR}/case.$1.$2 -res $1 -compset $compset -mach ${CCSM_MACH} > test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: create_newcase was successful" 
else
    echo "TCB_ccsm.sh: create_newcase failed, error from create_newcase= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

cd ${CAM_TESTDIR}/case.$1.$2
./xmlchange -file env_build.xml -id EXEROOT -val ${CAM_TESTDIR}/case.$1.$2 -silent
./xmlchange -file env_run.xml -id RUNDIR -val ${CAM_TESTDIR}/case.$1.$2/run -silent

# chemistry preprocessor
if [ $usrmech != $2 ]; then
   string1=`grep CAM_CONFIG_OPTS env_build.xml`
   string2=`echo $string1 | cut -d "=" -f 3`
   cfgstring=`echo $string2 | cut -d "\"" -f 2`
   ./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -val "$cfgstring -usr_mech_infile ${CAM_SCRIPTDIR}/config_files/$usrmech" 
fi

#
# Override CESM pes layouts on yellowstone
#
if [ ${CCSM_MACH} = 'yellowstone' ]; then
   ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 64
   ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 64
   ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 64
   ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 64
   ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 64
   ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 64
   ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 64
   ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 64

   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val 1
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val 1
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val 1
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val 1
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val 1
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val 1
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val 1
fi


./cesm_setup >> ${CAM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: CESM configure was successful" 
else
    echo "TCB_ccsm.sh: CESM configure failed, error from configure= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

cp ${CAM_SCRIPTDIR}/nl_files/user_nl_cam .

buildscript=`ls *.build`
./$buildscript >> ${CAM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: CESM build was successful" 
else
    echo "TCB_ccsm.sh: CESM build failed, error from build= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

cd ${blddir}
echo "TCB_ccsm.sh: CESM configure and build test passed"
echo "PASS" > TestStatus
if [ $CAM_RETAIN_FILES != "TRUE" ]; then
    echo "TCB_ccsm.sh: removing some unneeded files to save disc space" 
    rm -rf ${CAM_TESTDIR}/case.$1.$2/atm
    rm -rf ${CAM_TESTDIR}/case.$1.$2/glc
    rm -rf ${CAM_TESTDIR}/case.$1.$2/ice
    rm -rf ${CAM_TESTDIR}/case.$1.$2/lnd
    rm -rf ${CAM_TESTDIR}/case.$1.$2/ocn
    rm -rf ${CAM_TESTDIR}/case.$1.$2/rof
    rm -rf ${CAM_TESTDIR}/case.$1.$2/wav
    rm -rf ${CAM_TESTDIR}/case.$1.$2/cpl
    rm -rf ${CAM_TESTDIR}/case.$1.$2/mct
    rm -rf ${CAM_TESTDIR}/case.$1.$2/pio
    rm -rf ${CAM_TESTDIR}/case.$1.$2/cesm
    rm -rf ${CAM_TESTDIR}/case.$1.$2/csm_share
    rm -rf ${CAM_TESTDIR}/case.$1.$2/lib
fi

exit 0
