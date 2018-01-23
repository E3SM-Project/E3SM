#!/bin/tcsh 
#
#SBATCH -p ec
#SBATCH --job-name dcmip4
#SBATCH --account=FY150001
#SBATCH -N 4
#SBATCH --time=0:10:00
#XXSBATCH -N 54
#XXSBATCH --time=12:00:00
#PBS -l walltime=10:00
#PBS -l nodes=4
#PBS -q acme
#
# hydrostatic x1:  4 nodes, 3.3min
#

set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?PBS_NNODES} ) then   # redsky
    if ( $PBS_ENVIRONMENT == PBS_BATCH ) cd $PBS_O_WORKDIR     
    set NCPU = $PBS_NNODES
endif
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

# hydrostatic
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30

#set case  = h-x1
#set case  = nh-x1
set case  = nh-x1000

\rm -f input.nl
mkdir restart


# run 15 days, eulerian,  (no restart)
set namelist = ${case}-erun3.nl            
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

# run 15 days, write a restart file
set namelist = ${case}-erun1.nl            
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

# restart, run 2 timesteps with no viscosity:
set namelist = ${case}-erun2.nl            
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

exit 0



