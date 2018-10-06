#!/bin/bash
#
#SBATCH --job-name d16-1-preqx 
#SBATCH --account=FY150001
#SBATCH -N 30
#SBATCH --time=0:60:00
#SBATCH -p ec
#PBS -q acme
#PBS -l walltime=30:00
#PBS -l nodes=25
#
# 25 nodes, 30min sufficient for all 5 runs
# 12 nodes, 10min for r400 an r100
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
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30



function run { 
local NCPU=$1
echo "NCPU = $NCPU"
namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
date
mpirun -np $NCPU $EXEC < input.nl
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
#prefix=ne120  ; run $NCPU       # 72 nodes, 2h Anvil:  6s
#prefix=ne256  ; run $NCPU       # 72 nodes, 2h Anvil:  55s
                                 # 50 nodes ? 
                                 # 25 nodes ?   
#prefix=ne512  ; run $NCPU       # 72 nodes killed.  should work with 
#prefix=ne1024  ; run $NCPU       # 









