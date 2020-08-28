#!/bin/bash
#
##SBATCH --job-name d21-theta
##SNL:
##SBATCH --account=FY150001
##SBATCH -p ec
##anvil:
#SBATCH --account=condo
##SBATCH -p acme-small
#SBATCH -p acme-medium
#SBATCH -N 13
#SBATCH --time=5:00:00


export OMP_NUM_THREADS=1
# compute number of MPI tasks                                                                                               
if [ -n "$SLURM_NNODES" ]; then
  NNODES=$SLURM_NNODES
else
  NNODES=1
fi

#EXEC=~/runhomme/gw/test_execs/theta-l-nlev60/theta-l-nlev60
EXEC=../../../test_execs/theta-l-nlev60/theta-l-nlev60


if false; then
source ~/mods-nov25

# hydrostatic theta
namelist=namelist-h.nl
srun -K -c 1 -N $SLURM_NNODES  $EXEC < $namelist
#srun -K -c 1 -n ${1}  $EXEC < $namelist
source ~/load-nlc

ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_11.nc  movies/hydro_dcmip2012_test2_11.nc
\mv -f dcmip2012_test2_1_T_t10.pdf   hydro_test2_1_T_t10.pdf

fi


source ~/mods-nov25

# nonhydrostatic theta
#namelist=namelist-nh.nl
namelist=ne1024.nl
srun -K -c 1 -N $SLURM_NNODES  $EXEC < $namelist
#srun -K -c 1 -n ${1}  $EXEC < $namelist
source ~/load-nlc

ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_11.nc  movies/nonhydro_dcmip2012_test2_11.nc
\mv -f dcmip2012_test2_1_T_t10.pdf   nonhydro_test2_1_T_t10.pdf

cp *pdf ~/web/gw21

date
