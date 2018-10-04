#!/bin/tcsh 
#
#SBATCH -p ec
#SBATCH --job-name dcmip4
#SBATCH --account=FY150001
#SBATCH -N 12
#SBATCH --time=1:30:00
#XXSBATCH -N 20
#XXSBATCH --time=5:00:00
#PBS -l walltime=20:00
#PBS -l nodes=20
#PBS -q acme
#
# on anvil, < 10min to run all simuilations on 20 nodes
#
#
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

# hydrostatic
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = h-x1.nl           
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/hydro-X1-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/hydro-X1-dcmip2012_test41.nc"'
\mv -f zeta.pdf hydro-X1-zeta.pdf
\mv -f ps.pdf hydro-X1-ps.pdf

# nonhydro X1000
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = nh-x1000.nl       
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X1000-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X1000-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X1000-zeta.pdf
\mv -f ps.pdf nonhydro-X1000-ps.pdf

# nonhydro X100
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = nh-x100.nl        
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X100-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X100-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X100-zeta.pdf
\mv -f ps.pdf nonhydro-X100-ps.pdf


# nonhydro X10
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = nh-x10.nl        
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X10-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X10-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X10-zeta.pdf
\mv -f ps.pdf nonhydro-X10-ps.pdf

# nonhydro X1
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = nh-x1.nl         
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X1-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X1-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X1-zeta.pdf
\mv -f ps.pdf nonhydro-X1-ps.pdf


date



