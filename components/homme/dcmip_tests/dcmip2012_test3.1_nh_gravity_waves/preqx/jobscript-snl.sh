#!/bin/tcsh 
#
#   Jobscript for launching dcmip2012 test3-1 on SNL clusters
#
#SBATCH -J d31-preqx          # job name
#SBATCH -N 4                  # total number of mpi tasks requested
#SBATCH -p ec                 # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00           # run time (hh:mm:ss)
#SBATCH --account=FY150001

set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?SLURM_NNODES} ) then   # skybridge
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif
date

#############################################################################
# preqx (hydrostatic
#############################################################################
set EXEC = ../../../test_execs/preqx-nlev20-interp/preqx-nlev20-interp  
set namelist = ./namelist-default.nl

\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

ncl plot_omega.ncl
ncl plot_theta.ncl

\mv -f test31_omega.pdf                 preqx_test31_omega.pdf
\mv -f dcmip2012_test3_theta_diff.pdf   preqx_test3_theta_diff.pdf  
\mv -f movies/dcmip2012_test31.nc       movies/preqx_dcmip2012_test31.nc 

date

