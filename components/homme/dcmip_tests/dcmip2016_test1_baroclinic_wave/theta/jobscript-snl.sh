#!/bin/tcsh 
#
#SBATCH --job-name d16-1-theta
#SBATCH --account=FY150001
#SBATCH -N 12
#SBATCH --time=0:20:00
#SBATCH -p ec

set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif
date


set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
\cp -f namelist-r400-dry.nl input.nl
mpirun -np $NCPU $EXEC < input.nl

date

