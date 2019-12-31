#!/bin/tcsh 
#
#SBATCH --job-name dcmip4
#XXSBATCH -p ec
#XXSBATCH --account=FY150001
#SBATCH -p acme-medium
#SBATCH --account=condo           
#SBATCH -N 12
#SBATCH --time=1:30:00
#
# on anvil, < 10min to run all simuilations on 20 nodes
#
#
set OMP_NUM_THREADS = 1


# hydrostatic
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = h-x1.nl           
cp -f $namelist input.nl
srun -K -c 1 -N $SLURM_NNODES  $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/hydro-X1-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/hydro-X1-dcmip2012_test41.nc"'
\mv -f zeta.pdf hydro-X1-zeta.pdf
\mv -f ps.pdf hydro-X1-ps.pdf
\mv -f ps-8.pdf hydro-X1-ps-8.pdf

# nonhydro X1000
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = nh-x1000.nl       
cp -f $namelist input.nl
srun -K -c 1 -N $SLURM_NNODES $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X1000-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X1000-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X1000-zeta.pdf
\mv -f ps.pdf nonhydro-X1000-ps.pdf
\mv -f ps-8.pdf nonhydro-X1000-ps-8.pdf

# nonhydro X100
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = nh-x100.nl        
cp -f $namelist input.nl
srun -K -c 1 -N $SLURM_NNODES $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X100-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X100-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X100-zeta.pdf
\mv -f ps.pdf nonhydro-X100-ps.pdf
\mv -f ps-8.pdf nonhydro-X100-ps-8.pdf


# nonhydro X10
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = nh-x10.nl        
cp -f $namelist input.nl
srun -K -c 1 -N $SLURM_NNODES $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X10-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X10-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X10-zeta.pdf
\mv -f ps.pdf nonhydro-X10-ps.pdf
\mv -f ps-8.pdf nonhydro-X10-ps-8.pdf

# nonhydro X1
date
set EXEC = ../../../test_execs/theta-l-nlev30/theta-l-nlev30
set namelist = nh-x1.nl         
cp -f $namelist input.nl
srun -K -c 1 -N $SLURM_NNODES $EXEC < input.nl
ncl plot_ps.ncl  'fname="movies/nonhydro-X1-dcmip2012_test41.nc"'
ncl plot_zeta.ncl 'fname="movies/nonhydro-X1-dcmip2012_test41.nc"'
\mv -f zeta.pdf nonhydro-X1-zeta.pdf
\mv -f ps.pdf nonhydro-X1-ps.pdf
\mv -f ps-8.pdf nonhydro-X1-ps-8.pdf


date



