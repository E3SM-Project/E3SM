#!/bin/bash
# Script to test run_acme.template.csh
#   Absolute path to RUN_ACME_DIR is supplied by CMake configuration
#   Includes Batch-system specific method for checking queue
#      Supports:  pbs, slurm, none
#

#Case directories  are in HOME unless SCRATCH is specified
if [ -z "$SCRATCH" ]; then
  SCRATCH=${HOME}
fi

if [ ! -d "${SCRATCH}/ACME_simulations" ]; then
  mkdir ${SCRATCH}/ACME_simulations
fi

case_name="default.default.A_WCYCL1850_template.ne4np4_oQU240"
case_scratch="${SCRATCH}/ACME_simulations/${case_name}"
case_dir="${case_scratch}/case_scripts"
public_scratch="${SCRATCH}"
run_acme_log="${public_scratch}/run_acme.log"

#Verify the case is already deleted
rm -rf ${case_scratch}

# Change to ACME dir, where script must be run, using CMake for absolute path
cd @RUN_ACME_DIR@

# Run run_acme
echo "CTEST_FULL_OUTPUT" #This magic string stops CDash from truncating output
echo
echo "**********************************************"
echo "Running run_acme.template.csh :"
echo

./run_acme.template.csh > ${run_acme_log}
exitcode=$?
cat ${run_acme_log}
# Verify the script didn't fail
if [ $exitcode -ne 0 ]
  then exit $exitcode
fi

# cat CaseStatus so it is seen on CDash
echo
echo "**********************************************"
echo "cat ${case_dir}/CaseStatus :"
echo
cat "${case_dir}/CaseStatus"

# We want to wait for the last job to have been submitted to be complete
# Note that this assumes that only one job is submitted per case.submit;
# any more would require us to loop and wait for each of them
echo
echo "**********************************************"
echo "Wait for job to finish, if using batch system"
echo
jobid=`grep "Submitted job id is " ${run_acme_log} | cut -d' ' -f5 | tail -n1`

# CaseStatus doesn't report when the job has finished yet; let the queue tell us instead
# Query batch system name from case directory env_batch.xml
pushd ${case_dir}
batchsystem=$(./xmlquery BATCH_SYSTEM | sed "s/ //g" | sed "s/	//g" | sed "s/BATCH_SYSTEM://g")
popd

# Wait for job to finish using batch_system specific command
if [ "$batchsystem" == "none" ]; then
  echo "BATCH_SYSTEM = none ; No need to wait for batch job to finish"

elif [ "$batchsystem" == "pbs" ]; then
  echo "BATCH_SYSTEM = pbs ; Waiting for batch job to finish"
  while [ `qstat ${jobid} | wc -l` -eq 3 ]
    do sleep 5
  done

elif [ "$batchsystem" == "slurm" ]; then
  echo "BATCH_SYSTEM = slurm ; Waiting for batch job to finish"
  while [ `squeue -j ${jobid} | wc -l` -eq 2 ]
    do sleep 5
  done

else
  echo "ERROR: Unrecognized BATCH_SYSTEM = "$batchsystem
  echo "   Expected values are: none, pbs, slurm"
  echo "Not waiting for batch job to finish before checking for success"
fi
# 
# Paranoid delay to ensure we get the most up to date acme and coupler logs
sleep 5
sync

# Judging success by existence of specific string in cpl.log
#   The log files are typically gzipped, so use less to read them
find ${case_scratch}/run/ -iname "cpl.log.*" | xargs less | grep "SUCCESSFUL TERMINATION"
exitcode=$?
if [ $exitcode -ne 0 ]
  then exit $exitcode
fi
