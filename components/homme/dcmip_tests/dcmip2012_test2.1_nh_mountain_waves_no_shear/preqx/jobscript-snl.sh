#!/bin/tcsh 
#
#SBATCH --job-name d21
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

# NH model
set EXEC = ../../../test_execs/theta-nlev60/theta-nlev60    
set namelist = namelist-nh-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_11.nc  movies/nonhydro_dcmip2012_test2_11.nc
\mv -f dcmip2012_test2_1_T_t10.pdf   nonhydro_T_t10.pdf


# hydrostatic theta
set EXEC = ../../../test_execs/theta-nlev60/theta-nlev60    
set namelist = namelist-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_11.nc  movies/hydro_dcmip2012_test2_11.nc
\mv -f dcmip2012_test2_1_T_t10.pdf   hydro_T_t10.pdf


# hydrostatic preqx
set EXEC = ../../../test_execs/preqx-nlev60-interp/preqx-nlev60-interp        # set name of executable
set namelist = namelist-default.nl
\cp -f $namelist input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_11.nc  movies/preqx_dcmip2012_test2_11.nc
\mv -f dcmip2012_test2_1_T_t10.pdf   preqx_T_t10.pdf

date
