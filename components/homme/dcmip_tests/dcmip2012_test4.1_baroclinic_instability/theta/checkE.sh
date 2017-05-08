#!/bin/tcsh 
#
#SBATCH --job-name dcmip4
#SBATCH --account=FY150001
#SBATCH -N 4
#SBATCH --time=0:10:00
#XXSBATCH -N 54
#XXSBATCH --time=12:00:00
#SBATCH -p ec
#
# hydrostatic x1:  4 nodes, 3.3min
#

set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

# hydrostatic
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30

#set case  = h-x1
#set case  = nh-x100
set case  = nh-x1000

\rm -f input.nl
mkdir restart


# run 30 days, eulerian,  (no restart)
set namelist = ${case}-erun3.nl            
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

# run 30 days, write a restart file
set namelist = ${case}-erun1.nl            
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

# restart, run 2 timesteps with no viscosity:
set namelist = ${case}-erun2.nl            
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl




