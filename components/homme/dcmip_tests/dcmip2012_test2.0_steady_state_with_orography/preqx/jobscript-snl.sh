#!/bin/tcsh 
#
# nonhydrostatic model with explicit timestep: 25 nodes, 2h 20min
# hydrostatic: 300x faster.  4 nodes: 3min
#
#SBATCH --job-name d20
#SBATCH --account=FY150001
#SBATCH -N 25
#SBATCH --time=3:00:00
#XXSBATCH -N 4
#XXSBATCH --time=0:10:00
#SBATCH -p ec


set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?SLURM_NNODES} ) then   
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

# hydrostatic preqx
set EXEC = ../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
set namelist = namelist-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_z_lon.ncl
ncl test200-range.ncl
\mv -f dcmip2012_test2_0_u_t6.00.pdf preqx_test2_0_u_z.pdf
\mv -f movies/dcmip2012_test2_01.nc.pdf  preqx_test2_0_u.pdf
\mv -f movies/dcmip2012_test2_01.nc  movies/preqx_dcmip2012_test2_01.nc 



# hydrostatic theta
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30    
set namelist = namelist-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_z_lon.ncl
ncl test200-range.ncl
\mv -f dcmip2012_test2_0_u_t6.00.pdf hydro_test2_0_u_z.pdf
\mv -f movies/dcmip2012_test2_01.nc.pdf  hydro_test2_0_u.pdf
\mv -f movies/dcmip2012_test2_01.nc  movies/hydro_dcmip2012_test2_01.nc 



# NH model
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
set namelist = namelist-nh-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_z_lon.ncl
ncl test200-range.ncl
\mv -f dcmip2012_test2_0_u_t6.00.pdf nonhydro_test2_0_u_t6.00.pdf
\mv -f movies/dcmip2012_test2_01.nc.pdf  nonhydro_test2_0_u.pdf
\mv -f movies/dcmip2012_test2_01.nc  movies/nonhydro_dcmip2012_test2_01.nc 









