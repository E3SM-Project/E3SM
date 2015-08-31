#!/bin/sh
#==============================================================================
#
#  This script defines how to run CTest on the Argonne Leadership Computing
#  Facility systems (mira/cetus/vesta/cooley).
#
#  This assumes the CTest model name (e.g., "Nightly") is passed to it when
#  run.
#
#==============================================================================

# Get the CTest script directory
scrdir=$1

# Get the CTest model name
model=$2

# Write QSUB submission script with the test execution command
echo "#!/bin/sh" > runctest.pbs
echo "#PBS -q debug" >> runctest.pbs
echo "#PBS -l mppwidth=24" >> runctest.pbs
echo "#PBS -l walltime=00:20:00" >> runctest.pbs
echo "#PBS -V" >> runctest.pbs
echo "cd \$PBS_O_WORKDIR" >> runctest.pbs
echo "CTESTCMD=`which ctest`" >> runctest.pbs
echo "\$CTESTCMD -S ${scrdir}/CTestScript-Test.cmake,${model} -V" >> runctest.pbs

# Submit the job to the queue
jobid=`qsub runctest.pbs`

# Wait for the job to complete before exiting
while true; do
	status=`qstat $jobid`
	if [ "$status" == "" ]; then
		break
	else
		sleep 10
	fi
done
