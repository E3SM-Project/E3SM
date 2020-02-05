#!/bin/bash
#
#SBATCH --job-name d16-1-preqx 
#XXSBATCH --account=FY150001
#XXSBATCH -p ec
#SBATCH --account=condo
#SBATCH -p acme-medium
#SBATCH -N 25
#SBATCH --time=0:30:00
#
# 25 nodes, 30min sufficient for all 5 runs
# 12 nodes, 10min for r400 an r100
# 

export OMP_NUM_THREADS=1
export OMP_STACKSIZE=16M     #  Cori has 96GB per node. had to lower to 8M on 3K nodes
export MV2_ENABLE_AFFINITY=0
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

# hydrostatic preqx
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30



function run { 
local NCPU=$1
echo "NCPU = $NCPU"
namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
date
srun -K -c 1 -n $NCPU -N $SLURM_NNODES $EXEC < input.nl
date

ncl plot-baroclinicwave-init.ncl
ncl plot-lat-lon-TPLSPS.ncl 'var_choice=1'
ncl plot-lat-lon-TPLSPS.ncl 'var_choice=2'
ncl plot-lat-lon-TPLSPS.ncl 'var_choice=3'
\mv -f plot_baroclinicwave_init.pdf  ${prefix}_init.pdf
\mv -f preqx-test16-1latlonT850.pdf  ${prefix}_T850.pdf
\mv -f preqx-test16-1latlonPS.pdf  ${prefix}_PS.pdf
\mv -f preqx-test16-1latlonPRECL.pdf  ${prefix}_PRECL.pdf

\mv -f movies/dcmip2016_test11.nc    movies/${prefix}_dcmip2016_test11.nc
}

prefix=r400    ; run $(($NCPU>384?384:NCPU))

prefix=r100-dry; run $NCPU
prefix=r100-h  ; run $NCPU
prefix=r100    ; run $NCPU

prefix=r50    ; run $NCPU

# high res cases
#prefix=ne120  ; run $NCPU       
#prefix=ne256  ; run $NCPU       
#prefix=ne512  ; run $NCPU       
#prefix=ne1024  ; run $NCPU      

# timings on ANVIL
#ne120  72 nodes, 2h Anvil:  6s
#ne256  72 nodes, 2h Anvil:  55s
#ne512  100 nodes 2h Anvil:  611s (can run on 25 nodes)
#ne1024 100 nodes 2h Anvi:   4642s  (6min init)
#       100 nodes  270 timesteps:  1836s (6min init)
#



