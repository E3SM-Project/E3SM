#!/bin/bash
#
#SBATCH --job-name d22-theta
#SNL:
#SBATCH --account=FY150001
#SBATCH -p ec
#anvil:
#SBATCH --account=condo
#SBATCH -p acme-medium
#SBATCH -N 12
#SBATCH --time=0:20:00

export OMP_NUM_THREADS=1
# compute number of MPI tasks                                                                                               
if [ -n "$SLURM_NNODES" ]; then
  NNODES=$SLURM_NNODES
else
  NNODES=1
fi


EXEC=../../../test_execs/theta-l-nlev60/theta-l-nlev60

# hydrostatic theta
namelist=namelist-h.nl
\cp -f $namelist input.nl
srun -K -c 1 -N $SLURM_NNODES  $EXEC < input.nl
ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_21.nc  movies/hydro_dcmip2012_test2_21.nc
\mv -f dcmip2012_test2_2_T_t10.pdf   hydro_test2_2_T_t10.pdf

# nonhydrostatic theta
namelist=namelist-nh.nl
\cp -f $namelist input.nl
srun -K -c 1 -N $SLURM_NNODES   $EXEC < input.nl
ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_21.nc  movies/nonhydro_dcmip2012_test2_21.nc
\mv -f dcmip2012_test2_2_T_t10.pdf   nonhydro_test2_2_T_t10.pdf

date
