#!/bin/tcsh 
#
#   Jobscript for launching dcmip2012 test3-1 on SNL clusters
#
#SBATCH -J d31-theta          # job name
#SBATCH -N 4                  # total number of mpi tasks requested
#SBATCH -p ec                 # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00           # run time (hh:mm:ss)
#SBATCH --account=FY150001
#PBS -q acme
#PBS -l walltime=10:00                                                                                           
#PBS -l nodes=4    

set OMP_NUM_THREADS = 1
set NCPU = 40

if ( ${?PBS_ENVIRONMENT} ) then   # anvil
  set NCPU = $PBS_NNODES
  if ( $PBS_ENVIRONMENT == PBS_BATCH ) cd $PBS_O_WORKDIR     
endif
if ( ${?SLURM_NNODES} ) then   # skybridge
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

set EXEC = ../../../test_execs/theta-l-nlev20/theta-l-nlev20

#############################################################################
# theta (hydrostatic
#############################################################################
set namelist = ./namelist-h.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

ncl plot_omega.ncl
ncl plot_theta.ncl

\mv -f test31_omega.pdf                 hydro_test31_omega.pdf
\mv -f dcmip2012_test3_theta_diff.pdf   hydro_test3_theta_diff.pdf
\mv -f movies/dcmip2012_test31.nc       movies/hydro_dcmip2012_test31.nc

#############################################################################
# theta-nh
#############################################################################
set namelist = ./namelist-nh.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

ncl plot_omega.ncl
ncl plot_theta.ncl

\mv -f test31_omega.pdf                 nonhydro_test31_omega.pdf
\mv -f dcmip2012_test3_theta_diff.pdf   nonhydro_test3_theta_diff.pdf  
\mv -f movies/dcmip2012_test31.nc        movies/nonhydro_dcmip2012_test31.nc 
