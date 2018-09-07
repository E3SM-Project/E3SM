#!/bin/tcsh 
#
#SBATCH --job-name d16-1-preqx 
#SBATCH --account=FY150001
#SBATCH -N 12
#SBATCH --time=0:20:00
#SBATCH -p ec
#PBS -q acme
#PBS -l walltime=20:00
#PBS -l nodes=12


set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?PBS_ENVIRONMENT} ) then   # anvil
  set NCPU = $PBS_NNODES
  if ( $PBS_ENVIRONMENT == PBS_BATCH ) cd $PBS_O_WORKDIR     
endif
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif
date

# hydrostatic preqx
set EXEC = ../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
set prefix  = r100-dry
set namelist = namelist-$prefix.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl

ncl plot-baroclinicwave-init.ncl
ncl plot-lat-lon-TPLSPS.ncl
\mv -f plot_baroclinicwave_init.pdf  ${prefix}_init.pdf
\mv -f acme-test16-1latlonT850.pdf  ${prefix}_T850.pdf

\mv -f movies/dcmip2016_test11.nc    movies/${prefix}_dcmip2012_test11.nc


date
