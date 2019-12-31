#!/bin/bash -i
#
#SBATCH --job-name d16-2-theta-l
#SBATCH --account=FY150001
#SBATCH -p ec
#SBATCH --account=condo
#SBATCH -p acme-medium
#SBATCH -N 12
#SBATCH --time=0:20:00
#
#  r100 (ne30) 12 nodes, needs 5 min?
#  


OMP_NUM_THREADS=1
NCPU=8
if [ -n "$SLURM_NNODES" ]; then
   patt='([[:digit:]]+)'
   if [[ $SLURM_TASKS_PER_NODE =~ $patt ]]; then
     PER_NODE=${BASH_REMATCH[1]}
   else
     PER_NODE=16
   fi
   NCPU=$(( $SLURM_NNODES * $PER_NODE ))
fi

# theta model
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30


function run {
local NCPU=$1
echo "NCPU = $NCPU"
namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
date
srun -K -c 1 -n $NCPU -N $SLURM_NNODES  $EXEC < input.nl
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
\mv -f psmap.pdf ${prefix}_psmap.pdf
}


prefix=r400 ;  run $(($NCPU>384?384:NCPU))
prefix=r100 ;  run $NCPU
prefix=r50 ;  run $NCPU
prefix=r50-h ;  run $NCPU
