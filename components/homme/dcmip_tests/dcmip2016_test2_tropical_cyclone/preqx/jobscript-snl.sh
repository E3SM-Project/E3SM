#!/bin/bash -i
#
#SBATCH --job-name d16-2-preqx 
#SBATCH --account=FY150001
#SBATCH -N 12
#SBATCH --time=0:20:00
#SBATCH -p ec
#PBS -q acme
#PBS -l walltime=20:00
#PBS -l nodes=12
#
#  r100 (ne30) 12 nodes, needs 5 min?
#  


OMP_NUM_THREADS=1
NCPU=40
if [ -n "$PBS_ENVIRONMENT" ]; then
#  NCPU=$PBS_NNODES
  [ "$PBS_ENVIRONMENT" = "PBS_BATCH" ] && cd $PBS_O_WORKDIR 
  NCPU=$PBS_NNODES
fi
if [ -n "$SLURM_NNODES" ]; then
    NCPU=$SLURM_NNODES
    let NCPU*=16
    let NCPU/=$OMP_NUM_THREADS
fi

# hydrostatic preqx
EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp  


function run {
local NCPU=$1
echo "NCPU = $NCPU"
namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
date
mpirun -np $NCPU $EXEC < input.nl
date

ncl plot-tropical-cyclone-init.ncl  # u,t,th,q,pnh,geo,ps, time=0
ncl plot-horiz-crossx.ncl     # contour plot, time=10d, U,V,T,ps,precl,Q,geo
ncl plot-intensity-trace.ncl
ncl plot-horiz-ps.ncl

# save output
\mv -f movies/dcmip2016_test21.nc   movies/${prefix}_dcmip2016_test21.nc

\mv -f init.pdf ${prefix}_init.pdf
\mv -f x-sections.pdf ${prefix}_x-sections.pdf
\mv -f wind.pdf ${prefix}_wind.pdf
\mv -f ps.pdf ${prefix}_ps.pdf
\mv -f psmap.pdf ${prefix}_psmap.pdf
}


prefix=r400 ;  run $(($NCPU>384?384:NCPU))
prefix=r100 ;  run $NCPU
prefix=r50 ;  run $NCPU
