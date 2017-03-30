#!/bin/tcsh 
#
# hydrostatic: 4 nodes: 3min
#
#SBATCH --job-name d20-preqx
#SBATCH --account=FY150001
#SBATCH -N 4
#SBATCH --time=0:10:00
#SBATCH -p ec

set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?SLURM_NNODES} ) then   
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

# hydrostatic preqx
set EXEC = ../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp  # set name of executable
set namelist = namelist-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_z_lon.ncl
ncl test200-range.ncl







