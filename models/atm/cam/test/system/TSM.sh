#!/bin/sh 
#

echo $1 $2 $3 $4

if [ $# -ne 4 ]; then
    echo "TSM.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSM.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSM.sh: smoke test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSM.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSM.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CAM_TESTDIR}/${test_name} ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=${CAM_SCRIPTDIR}/../../bld
rundir=${CAM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -rf ${rundir}
fi
mkdir -p ${rundir}
if [ $? -ne 0 ]; then
    echo "TSM.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSM.sh: calling TCB.sh to prepare cam executable" 
${CAM_SCRIPTDIR}/TCB.sh $1 $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSM.sh: error from TCB.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ $4 = "build_only" ]; then
  exit 0
fi 

nl_file=${2%+*}
use_case=${2#*+}
cfg_file=${1%+*}

if [ ! -f ${CAM_SCRIPTDIR}/nl_files/$nl_file ]; then
    echo "TSM.sh: namelist options file ${CAM_SCRIPTDIR}/nl_files/$nl_file not found" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

if [ ${use_case} != $2 ]; then
    use_case_string="-use_case ${use_case}"
else
    use_case_string=""
fi

stop_flag=${3##*[0-9]}
run_length=${3%%[sdm]}
if [ ${stop_flag} = $3 ] || [ ${run_length} = $3 ]; then
    echo "TSM.sh: error processing input argument for run length= $3"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 99
fi

case $stop_flag in
    s )  stop_option="nsteps";;

    d )  stop_option="ndays";;

    m )  stop_option="nmonths";;
esac

echo "TSM.sh calling CAM_decomp.sh to build decomposition string" 
${CAM_SCRIPTDIR}/CAM_decomp.sh $cfg_file $nl_file 1
rc=$?
if [ $rc -eq 0 ] && [ -f cam_decomp_string.txt ]; then
    read decomp_str < cam_decomp_string.txt
    echo "TSM.sh: cam decomp string: $decomp_str"  
    rm cam_decomp_string.txt
else
    echo "TSM.sh: error building decomp string; error from CAM_decomp.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

cp ${CAM_SCRIPTDIR}/nl_files/$nl_file ${CAM_TESTDIR}/${test_name}/.
perl -pi -e 's/\$CSMDATA/$ENV{CSMDATA}/'  ${CAM_TESTDIR}/${test_name}/$nl_file

#
# Turn on all outputs, when physics is cam3, cam4, or cam5
#

if [[ "$1" == c3 || "$1" == c4 || "$1" == c5 ]]; then
  history_output="history_aerosol=.true.  history_aero_optics=.true.  history_eddy=.true.  history_budget=.true."
fi

echo "TSM.sh: running cam; output in ${CAM_TESTDIR}/${test_name}/test.log" 
echo "TSM.sh: call to build-namelist:" 
echo "        env OMP_NUM_THREADS=${CAM_THREADS} ${cfgdir}/build-namelist -test -runtype startup \
    -ignore_ic_date $use_case_string \
    -config ${CAM_TESTDIR}/TCB.$1/config_cache.xml \
    -config_cice ${CAM_TESTDIR}/TCB.$1/config_cache_cice.xml \
    -infile ${CAM_TESTDIR}/${test_name}/$nl_file \
    -namelist \"&seq_timemgr_inparm stop_n=$run_length stop_option=\'$stop_option\' $decomp_str /\""  

env OMP_NUM_THREADS=${CAM_THREADS} ${cfgdir}/build-namelist -test -runtype startup -cice_nl "&domain_nml distribution_type='roundrobin' /" \
    -ignore_ic_date $use_case_string \
    -config ${CAM_TESTDIR}/TCB.$1/config_cache.xml \
    -config_cice ${CAM_TESTDIR}/TCB.$1/config_cache_cice.xml \
    -infile ${CAM_TESTDIR}/${test_name}/$nl_file \
    -namelist "&seq_timemgr_inparm stop_n=$run_length stop_option='$stop_option' $decomp_str $history_output/"  > test.log 2>&1
rc=$?
cat >ocn_modelio.nml <<EOF
&modelio
 diro = '$PWD'
 logfile = 'ocn.log'
/
&pio_inparm
 pio_numiotasks = -99
 pio_root = -99
 pio_stride = -99
 pio_typename = 'nothing'
/
EOF

if [ $rc -eq 0 ]; then
    echo "TSM.sh: cam build-namelist was successful" 
    cat drv_in
    cat atm_in
    cat lnd_in
    cat ocn_in
    cat ice_in
    cat docn_in
    cat docn_ocn_in
    cat drv_flds_in
    cat docn.stream.txt
    cat ocn_modelio.nml
else
    echo "TSM.sh: error building namelist, error from build-namelist= $rc" 
    echo "TSM.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

echo "TSM.sh calling CAM_runcmnd.sh to build run command" 
${CAM_SCRIPTDIR}/CAM_runcmnd.sh $cfg_file 1
rc=$?
if [ $rc -eq 0 ] && [ -f cam_run_command.txt ]; then
    read cmnd < cam_run_command.txt
    echo "TSM.sh: cam run command:" 
    echo "        $cmnd ${CAM_TESTDIR}/TCB.$1/cam"  
    rm cam_run_command.txt
else
    echo "TSM.sh: error building run command; error from CAM_runcmnd.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

${cmnd} ${CAM_TESTDIR}/TCB.$1/cam >> test.log 2>&1
rc=$?
if [ $rc -eq 0 ] && grep -c "END OF MODEL RUN" test.log > /dev/null; then
    echo "TSM.sh: smoke test passed" 
    echo "PASS" > TestStatus
    if [ $CAM_RETAIN_FILES != "TRUE" ]; then
        echo "TSM.sh: removing some unneeded files to save disc space" 
        if [ -f *.clm*.i.* ]; then
            rm *.clm*.i.*
	fi
    fi
else
    echo "TSM.sh: error running cam, error= $rc" 
    echo "TSM.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

exit 0
