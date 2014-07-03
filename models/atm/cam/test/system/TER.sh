#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TER.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TER.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TER.sh: exact restart test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TER.sh: exact restart test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TER.sh: this exact restart test failed under job ${prev_jobid} - moving those results to "
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
    echo "TER.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

initial_length=${3%+*}
restart_string=${3#*+}
if [ ${initial_length} = $3 ] || [ ${restart_string} = $3 ]; then
    echo "TER.sh: error processing input argument for run lengths"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

stop_flag=${restart_string##*[0-9]}
restart_length=${restart_string%%[sdm]}
if [ ${stop_flag} = ${restart_string} ] || [ ${restart_length} = ${restart_string} ]; then
    echo "TER.sh: error processing input argument for run length= $3"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

case $stop_flag in
    s )  stop_option="nsteps";;

    d )  stop_option="ndays";;

    m )  stop_option="nmonths";;
esac

full_length=`expr $initial_length + $restart_length`

echo "TER.sh: calling TSM.sh for smoke test of full length ${full_length}${stop_flag}" 
${CAM_SCRIPTDIR}/TSM.sh $1 $2 "${full_length}${stop_flag}" $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

if [ $4 = "build_only" ]; then
  exit 0
fi

echo "TER.sh: calling TSM.sh for smoke test of initial length ${initial_length}${stop_flag}" 
${CAM_SCRIPTDIR}/TSM.sh $1 $2 "${initial_length}${stop_flag}" $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi


cp ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*cam* ${rundir}/.
cp ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*clm* ${rundir}/.
cp ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*drv* ${rundir}/.
cp ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/rpointer* ${rundir}/.

use_case=${2#*+}
nl_file=${2%+*}

if grep -ic endofrun  ${CAM_SCRIPTDIR}/nl_files/$nl_file > /dev/null; then
  # REMOVE the inithist file created at the end of the initial_length run as it
  # will not be produced by the full_length run and will cause the test to fail.
  echo "TER.sh: Removing inithist written by initial_length run"
  rm ${rundir}/*cam*.i.*nc 
fi


if [ ${use_case} != $2 ]; then
    use_case_string="-use_case ${use_case}"
else
    use_case_string=""
fi

## modify the # of tasks/threads for restart if not testing fv decomposition and not a cosp run
if grep -ic npr_yz ${CAM_SCRIPTDIR}/nl_files/$nl_file > /dev/null || grep -ic cosp ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
    echo "TER.sh: will not modify tasks/threads for restart"
    factor=1
else
    echo "TER.sh: will modify tasks/threads for restart to tasks=${CAM_RESTART_TASKS} and threads=${CAM_RESTART_THREADS}"
    export CAM_TASKS=${CAM_RESTART_TASKS}
    export CAM_THREADS=${CAM_RESTART_THREADS}
    factor=2
fi

echo "TER.sh calling CAM_decomp.sh to build decomposition string"
${CAM_SCRIPTDIR}/CAM_decomp.sh $1 $nl_file $factor
rc=$?
if [ $rc -eq 0 ] && [ -f cam_decomp_string.txt ]; then
    read decomp_str < cam_decomp_string.txt
    echo "TER.sh: cam decomp string: $decomp_str"
    rm cam_decomp_string.txt
else
    echo "TER.sh: error building decomp string; error from CAM_decomp.sh= $rc"
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

echo "TER.sh: restarting sequential ccsm; output in ${CAM_TESTDIR}/${test_name}/test.log" 
echo "TER.sh: call to build-namelist:"
echo "        env OMP_NUM_THREADS=${CAM_THREADS} ${cfgdir}/build-namelist -test -runtype continue \
    -config ${CAM_TESTDIR}/TCB.$1/config_cache.xml \
    -config_cice ${CAM_TESTDIR}/TCB.$1/config_cache_cice.xml \
    -infile ${CAM_SCRIPTDIR}/nl_files/$nl_file \
    -infile ${CAM_TESTDIR}/${test_name}/$nl_file \
    -namelist \"&seq_timemgr_inparm stop_n=${restart_length} stop_option=\'$stop_option\' $decomp_str $history_output /\""  

env OMP_NUM_THREADS=${CAM_THREADS} ${cfgdir}/build-namelist -test -runtype continue -cice_nl "&domain_nml distribution_type='roundrobin' /" \
    -config ${CAM_TESTDIR}/TCB.$1/config_cache.xml \
    -config_cice ${CAM_TESTDIR}/TCB.$1/config_cache_cice.xml \
    -ignore_ic_date $use_case_string \
    -infile ${CAM_TESTDIR}/${test_name}/$nl_file \
    -namelist "&seq_timemgr_inparm stop_n=${restart_length} stop_option='$stop_option' $decomp_str $history_output / &cam_inparm /" > test.log 2>&1
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
    echo "TER.sh: cam build-namelist was successful"
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
    echo "TER.sh: error building namelist, error from build-namelist= $rc"
    echo "TER.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

echo "TER.sh calling CAM_runcmnd.sh to build run command"
${CAM_SCRIPTDIR}/CAM_runcmnd.sh $1 $factor
rc=$?
if [ $rc -eq 0 ] && [ -f cam_run_command.txt ]; then
    read cmnd < cam_run_command.txt
    echo "TER.sh: cam run command:"
    echo "        $cmnd ${CAM_TESTDIR}/TCB.$1/cam"
    rm cam_run_command.txt
else
    echo "TER.sh: error building run command; error from CAM_runcmnd.sh= $rc"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 9
fi

${cmnd} ${CAM_TESTDIR}/TCB.$1/cam >> test.log 2>&1
rc=$?
if [ $rc -eq 0 ] && grep -c "END OF MODEL RUN" test.log > /dev/null; then
    echo "TER.sh: restart of sequential ccsm completed successfully"
else
    echo "TER.sh: error on restart run of sequential ccsm, error= $rc" 
    echo "TER.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10
fi

echo "TER.sh: starting b4b comparisons " 
files_to_compare=`ls *.cam*.h*.nc *.cam*.i*.nc`
if [ -z "${files_to_compare}" ]; then
    echo "TER.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 11
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    ${CAM_SCRIPTDIR}/CAM_compare.sh \
	${compare_file} \
	${CAM_TESTDIR}/TSM.$1.$2.${full_length}${stop_flag}/${compare_file}
    rc=$?
    mv cprnc.out cprnc.${compare_file}.out
    if [ $rc -eq 0 ]; then
        echo "TER.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
	echo "TER.sh: error from CAM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	all_comparisons_good="FALSE"
    fi
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TER.sh: exact restart test passed" 
    echo "PASS" > TestStatus
    if [ $CAM_RETAIN_FILES != "TRUE" ]; then
        echo "TER.sh: removing some unneeded files to save disc space" 
        rm *.nc
        rm *.r*
    fi
else
    echo "TER.sh: at least one file comparison did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 12
fi

exit 0
