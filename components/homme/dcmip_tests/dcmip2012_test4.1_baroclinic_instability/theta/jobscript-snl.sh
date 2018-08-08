#!/bin/tcsh 
#
#SBATCH -p ec
#SBATCH --job-name dcmip4
#SBATCH --account=FY150001
#SBATCH -N 12
#SBATCH --time=1:30:00
#XXSBATCH -N 20
#XXSBATCH --time=5:00:00
#PBS -l walltime=60:00
#PBS -l nodes=20
#PBS -q acme

#
#  nonhydro x1:  54 nodes, 7.4h        KG5 dt=.5
#           x1:  20 nodes, 2.5h        ars232  dt=120
#           x10  54 nodes, 115min      KG5 dt=.2
#           x10  20 nodes, 2.5h        ars232 dt=12  (IN)
#           x100  12 nodes  51min
#           x1000 12 nodes 5.9min
# hydrostatic x1:   12 nodes, 3.3min

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
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
set namelist = h-x1.nl            #  4320 timesteps, 12 nodes, 3.3min
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/hydro-X1-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/hydro-X1-dcmip2012_test41.nc"'
\mv -f zeta.pdf hydro-X1-zeta.pdf
\mv -f ps.pdf hydro-X1-ps.pdf

# nonhydro X1000
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
set namelist = nh-x1000.nl        #  6480 timesteps, 12 nodes - 5.9min
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X1000-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X1000-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X1000-zeta.pdf
\mv -f ps.pdf nonhydro-X1000-ps.pdf

# nonhydro X100
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
set namelist = nh-x100.nl        
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X100-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X100-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X100-zeta.pdf
\mv -f ps.pdf nonhydro-X100-ps.pdf


# nonhydro X10
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
set namelist = nh-x10.nl         
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X10-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X10-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X10-zeta.pdf
\mv -f ps.pdf nonhydro-X10-ps.pdf

# nonhydro X1
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
set namelist = nh-x1.nl          
cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X1-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X1-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X1-zeta.pdf
\mv -f ps.pdf nonhydro-X1-ps.pdf





