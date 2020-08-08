#!/bin/bash 
#
# hydrostatic: 4 nodes: 3min
# NH:  10 nodes ?  
#
#SBATCH --job-name d20-theta
#XXSBATCH -p short,batch
#XXSBATCH --account=FY150001
#SBATCH -p acme-small
#SBATCH --account=condo
#SBATCH -N 4
#SBATCH --time=0:60:00


export OMP_NUM_THREADS=1
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev128    

source ~/mods-nov25 


# hydrostatic theta
namelist=namelist-h.nl
srun -K -c 1 -N $SLURM_NNODES  $EXEC < $namelist
source ~/load-nlc
ncl plot_z_lon.ncl
ncl test200-range.ncl
\mv -f dcmip2012_test2_0_u_t6.00.pdf hydro_test2_0_u_z.pdf
\mv -f movies/dcmip2012_test2_01.nc.pdf  hydro_test2_0_u.pdf
\mv -f movies/dcmip2012_test2_01.nc  movies/hydro_dcmip2012_test2_01.nc 

source ~/mods-nov25 

# nonhydrostatic theta
namelist=namelist-nh.nl
srun -K -c 1 -N $SLURM_NNODES  $EXEC < $namelist
source ~/load-nlc
ncl plot_z_lon.ncl
ncl test200-range.ncl
\mv -f dcmip2012_test2_0_u_t6.00.pdf nonhydro_test2_0_u_t6.00.pdf
\mv -f movies/dcmip2012_test2_01.nc.pdf  nonhydro_test2_0_u.pdf
\mv -f movies/dcmip2012_test2_01.nc  movies/nonhydro_dcmip2012_test2_01.nc 

cp *pdf ~/web/gw20



