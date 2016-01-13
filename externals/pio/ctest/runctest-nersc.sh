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
echo "#!/bin/sh" > runctest.slurm
echo "#SBATCH --partition debug" >> runctest.slurm
echo "#SBATCH --nodes=1" >> runctest.slurm
echo "#SBATCH --ntasks-per-node=32" >> runctest.slurm
echo "#SBATCH --time=00:15:00" >> runctest.slurm

echo "#SBATCH --export PIO_DASHBOARD_SITE,PIO_DASHBOARD_BUILD_NAME,PIO_DASHBOARD_SOURCE_DIR,PIO_DASHBOARD_BINARY_DIR" >> runctest.slurm
#echo "cd \$PBS_O_WORKDIR" >> runctest.pbs
echo "CTEST_CMD=`which ctest`" >> runctest.slurm
echo "\$CTEST_CMD -S ${scrdir}/CTestScript-Test.cmake,${model} -V" >> runctest.slurm
chmod +x runctest.slurm
# Submit the job to the queue
#jobid=`sbatch runctest.slurm| egrep -o -e "\b[0-9]+$"`
salloc -N 1 ./runctest.slurm
# Wait for the job to complete before exiting
#while true; do
#	status=`squeue -j $jobid`
#	if [ "$status" == "" ]; then
#		break
#	else
#		sleep 10
#	fi
#done
