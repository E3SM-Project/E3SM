#!/bin/tcsh 
#
#SBATCH --job-name d16-3-preqx 
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

# hydrostatic preqx
set EXEC = ../../../test_execs/preqx-nlev40-interp/preqx-nlev40-interp        # set name of executable
set namelist = namelist-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

date
