#!/bin/tcsh 
#
#SBATCH --job-name d22
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

set EXEC = ../../../test_execs/preqx-nlev60-interp/preqx-nlev60-interp      
set namelist = namelist-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_21.nc movies/preqx_dcmip2012_test2_21.nc
\mv -f dcmip2012_test2_2_T_t10.pdf preqx_T_t10.pdf

date
