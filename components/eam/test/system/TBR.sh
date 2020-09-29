#!/bin/sh 
#

if [ $# -ne 4 ]; then
    echo "TBR.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TBR.$1.$2.$3

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBR.sh: branch test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}"
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TBR.sh: branch test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TBR.sh: this branch test failed under job ${prev_jobid} - moving those results to "
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
    echo "TBR.sh: error, unable to create work subdirectory" 
    exit 3
fi

cd ${rundir}

initial_length=${3%+*}
branch_string=${3#*+}
if [ ${initial_length} = $3 ] || [ ${branch_string} = $3 ]; then
    echo "TBR.sh: error processing input argument for run lengths"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

stop_flag=${branch_string##*[0-9]}
branch_length=${branch_string%%[sdm]}
if [ ${stop_flag} = ${branch_string} ] || [ ${branch_length} = ${branch_string} ]; then
    echo "TBR.sh: error processing input argument for run length= $3"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

case $stop_flag in
    s )  stop_option="nsteps";;

    d )  stop_option="ndays";;

    m )  stop_option="nmonths";;
esac

#full_length=`expr $initial_length + $branch_length`
full_length=`expr $initial_length`

echo "TBR.sh: calling TSM.sh for smoke test of full length ${full_length}${stop_flag}"
${CAM_SCRIPTDIR}/TSM.sh $1 $2 "${full_length}${stop_flag}" $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBR.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

if [ $4 = "build_only" ]; then
  exit 0
fi

echo "TBR.sh: calling TSM.sh for smoke test of initial length ${initial_length}${stop_flag}"
${CAM_SCRIPTDIR}/TSM.sh $1 $2 "${initial_length}${stop_flag}" $4
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBR.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

cp ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*cam* ${rundir}/.

nl_file=${2%+*}
use_case=${2#*+}

if [ ${use_case} != $2 ]; then
    use_case_string="-use_case ${use_case}"
else
    use_case_string=""
fi

## branch from an older restart file if more than one available 
## note i need to account for clm generating additional netcdf restart files 
## note that there may be no clm file in certain modes (adiabatic for example) 
## and cam will ignore the namelist vars for clm
## the following commands now do their ls in the directory the files were 
## copied from because the above cp's sometimes wouldn't finish in time on bluesky
master_cam_restart=`ls -1rt ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*.cam*.r.*[0-9].nc \
    | tail -2 | head -1`
master_clm_restart=`ls -1rt ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*.clm*.r.* \
    | tail -3 | head -1`
master_cpl_restart=`ls -1rt ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*.cpl*.r.* \
    | tail -2 | head -1`
if [ -f ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/docn_in ]; then
   ocn_inparm=docn_in
   ocn_branch_nlparm=restfilm
   master_ocn_restart=`ls -1rt ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*.docn*.rs1.* \
       | tail -2 | head -1`
else
   ocn_inparm=dom_inparm
   ocn_branch_nlparm=dom_branch_file
   master_ocn_restart=`ls -1rt ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*.camdom*.r.* \
       | tail -2 | head -1`
fi
master_cice_path=`ls -1rt ${CAM_TESTDIR}/TSM.$1.$2.${initial_length}${stop_flag}/*.cice.r.[0-9]* \
    | tail -2 | head -1`
# cice_ic requires just the filename, not the absolute path
master_cice_restart=${master_cice_path##/*/}

## modify the # of tasks/threads for restart if not testing fv decomposition
if grep -ic npr_yz ${CAM_SCRIPTDIR}/nl_files/$nl_file > /dev/null; then
    echo "TBR.sh: configured for fv2d, will not modify tasks/threads for branch"
else
    echo "TBR.sh: will modify tasks/threads for branch to tasks=${CAM_RESTART_TASKS} and threads=${CAM_RESTART_THREADS}"
    export CAM_TASKS=${CAM_RESTART_TASKS}
    export CAM_THREADS=${CAM_RESTART_THREADS}
fi

echo "TBR.sh calling CAM_decomp.sh to build decomposition string"
${CAM_SCRIPTDIR}/CAM_decomp.sh $1 $nl_file 2
rc=$?
if [ $rc -eq 0 ] && [ -f cam_decomp_string.txt ]; then
    read decomp_str < cam_decomp_string.txt
    echo "TBR.sh: cam decomp string: $decomp_str"
    rm cam_decomp_string.txt
else
    echo "TBR.sh: error building decomp string; error from CAM_decomp.sh= $rc"
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

echo "TBR.sh: branching cam; output in ${CAM_TESTDIR}/${test_name}/test.log"
echo "TBR.sh: call to build-namelist:"
echo "        env OMP_NUM_THREADS=${CAM_THREADS} ${cfgdir}/build-namelist -test -runtype branch  \
    -config ${CAM_TESTDIR}/TCB.$1/config_cache.xml \
    -config_cice ${CAM_TESTDIR}/TCB.$1/config_cache_cice.xml \
    -infile ${CAM_TESTDIR}/${test_name}/$nl_file \
    -ignore_ic_date $use_case_string \
    -cice_nl \"&ice ice_ic=\'${master_cice_restart}\' /\" \
    -namelist \"&seq_timemgr_inparm stop_n=${branch_length} stop_option=\'$stop_option\' $decomp_str $history_output / \
    &seq_infodata_inparm brnch_retain_casename=.true. restart_file=\'${master_cpl_restart}\' / \
    &cam_inparm cam_branch_file=\'${master_cam_restart}\' / \
    &${ocn_inparm} ${ocn_branch_nlparm}='${master_ocn_restart} ' / \
    &clm_inparm nrevsn=\'${master_clm_restart}\' /\""

env OMP_NUM_THREADS=${CAM_THREADS} ${cfgdir}/build-namelist -test -runtype branch -cice_nl "&domain_nml distribution_type='roundrobin' /" \
    -config ${CAM_TESTDIR}/TCB.$1/config_cache.xml \
    -config_cice ${CAM_TESTDIR}/TCB.$1/config_cache_cice.xml \
    -ignore_ic_date $use_case_string \
    -infile ${CAM_TESTDIR}/${test_name}/$nl_file \
    -cice_nl "&ice ice_ic='${master_cice_restart}' /" \
    -namelist "&seq_timemgr_inparm stop_n=${branch_length} stop_option='$stop_option' $decomp_str $history_output / \
    &seq_infodata_inparm brnch_retain_casename=.true. \
    restart_file='${master_cpl_restart}' / \
    &cam_inparm cam_branch_file='${master_cam_restart}' / \
    &${ocn_inparm} ${ocn_branch_nlparm}='${master_ocn_restart} ' / \
    &clm_inparm nrevsn='${master_clm_restart}' /" \
     > test.log 2>&1
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
    echo "TBR.sh: cam build-namelist was successful"
    cat drv_in
    cat atm_in
    cat lnd_in
    cat docn_in
    cat ocn_in
    cat ice_in
    cat docn_in
    cat docn_ocn_in
    cat drv_flds_in
    cat docn.stream.txt
    cat ocn_modelio.nml
else
    echo "TBR.sh: error building namelist, error from build-namelist= $rc"
    echo "TBR.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

echo "TBR.sh calling CAM_runcmnd.sh to build run command"
${CAM_SCRIPTDIR}/CAM_runcmnd.sh $1 2
rc=$?
if [ $rc -eq 0 ] && [ -f cam_run_command.txt ]; then
    read cmnd < cam_run_command.txt
    echo "TBR.sh: cam run command:"
    echo "        $cmnd ${CAM_TESTDIR}/TCB.$1/cam"
    rm cam_run_command.txt
else
    echo "TBR.sh: error building run command; error from CAM_runcmnd.sh= $rc"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 9
fi

${cmnd} ${CAM_TESTDIR}/TCB.$1/cam >> test.log 2>&1
rc=$?
if [ $rc -eq 0 ] && grep -c "END OF MODEL RUN" test.log > /dev/null; then
    echo "TBR.sh: branch of cam completed successfully"
else
    echo "TBR.sh: error on branch run of cam, error= $rc"
    echo "TBR.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10
fi

echo "TBR.sh: starting b4b comparisons " 
files_to_compare=`ls *.cam*.h*.nc *.cam*.i*.nc`
if [ -z "${files_to_compare}" ]; then
    echo "TBR.sh: error locating files to compare"
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
        echo "TBR.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
	echo "TBR.sh: error from CAM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	all_comparisons_good="FALSE"
    fi
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TBR.sh: branch test passed" 
    echo "PASS" > TestStatus
    if [ $CAM_RETAIN_FILES != "TRUE" ]; then
        echo "TBR.sh: removing some unneeded files to save disc space" 
        rm *.nc
        rm *.r*
    fi
else
    echo "TBR.sh: at least one file comparison did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 12
fi

exit 0
