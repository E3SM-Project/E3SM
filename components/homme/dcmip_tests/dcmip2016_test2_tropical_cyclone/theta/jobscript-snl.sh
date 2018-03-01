#!/bin/csh 
#
#   Jobscript for launching dcmip2016 on SNL systems
#
#
#SBATCH --job-name d16-3-theta
#SBATCH --account=FY150001
#SBATCH -N 36
#SBATCH --time=0:30:00
#SBATCH -p ec

set NCPU = 40
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif



set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
#\cp ./namelist-r400.nl input.nl
\cp ./namelist-r400-nh.nl input.nl
mpirun -n $NCPU   $EXEC < input.nl                          
