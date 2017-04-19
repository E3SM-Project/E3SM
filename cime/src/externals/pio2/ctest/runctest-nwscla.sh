#!/bin/sh
#==============================================================================
#
#  This script defines how to run CTest on the NCAR Wyoming Supercomputing 
#  Center systems (yellowstone/caldera/geyser).
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
echo "#!/bin/sh" > runctest.sh
echo "#PBS -l walltime=01:00:00" >> runctest.sh
echo "#PBS -l select=1:ncpus=8:mpiprocs=8" >> runctest.sh
echo "#PBS -A SCSG0002" >> runctest.sh
echo "export PIO_DASHBOARD_SITE=nwscla-${HOSTNAME}" >> runctest.sh
echo "CTESTCMD=`which ctest`" >> runctest.sh
echo "\$CTESTCMD -S ${scrdir}/CTestScript-Test.cmake,${model} -V" >> runctest.sh

# Make the QSUB script executable
chmod +x runctest.sh

# Submit the job to the queue
jobid=`qsub runctest.sh`

# Wait for the job to complete before exiting
while true; do
	status=`qstat $jobid`
	echo $status
	if [ "$status" == "" ]; then
		break
	else
		sleep 10
	fi
done

exit 0

