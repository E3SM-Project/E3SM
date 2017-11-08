#!/bin/tcsh 
#
#SBATCH --job-name dcmip4
#SBATCH --account=FY150001
#SBATCH -N 12
#SBATCH --time=0:05:00
#SBATCH -p ec
#
# hydrostatic x1:   12 nodes, 3.3min
#

set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

#  preqx
set EXEC = ../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
set namelist = h-x1.nl      # 4320 timesteps,   12 nodes, 2.2min
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl ps.ncl   'fname="movies/preqx-X1-dcmip2012_test41.nc"'
ncl zeta.ncl 'fname="movies/preqx-X1-dcmip2012_test41.nc"'
\mv -f zeta.pdf preqx-X1-zeta.pdf
\mv -f ps.pdf preqx-X1-ps.pdf



