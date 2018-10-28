#!/bin/bash
#
# Anvil, 25 nodes, all 4 resolutions run in 25min
# need at least 11 nodes, since r400 is hardcoded as mpirun -np 384
#
#SBATCH --job-name d16-3-theta
#SBATCH --account=FY150001
#SBATCH -N 36
#SBATCH --time=0:30:00
#SBATCH -p ec
#PBS -q acme
#PBS -l walltime=30:00
#PBS -l nodes=25    

OMP_NUM_THREADS=1
NCPU=640
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

EXEC=../../../test_execs/theta-l-nlev40/theta-l-nlev40



function run { 
local NCPU=$1
echo "NCPU = $NCPU"
namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
date
mpirun -np $NCPU $EXEC < input.nl
date

ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

\mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_${prefix}.nc
\mv HommeTime                  HommeTime_${prefix}
\mv max_w.pdf                  max_w_${prefix}.pdf
\mv max_precip.pdf             max_precip_${prefix}.pdf
\mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_${prefix}.pdf
\mv 5km_xsec.pdf               5km_xsec_${prefix}.pdf
\mv measurement_wmax.txt       measurement_wmax_${prefix}.txt
\mv measurement_time.txt       measurement_time_${prefix}.txt
\mv measurement_prect_rate.txt measurement_prect_rate_${prefix}.txt

}

prefix=r400             ; run $(($NCPU>384?384:NCPU))
prefix=r200             ; run $NCPU
prefix=r100             ; run $NCPU
prefix=r50              ; run $NCPU

prefix=explicit-r100    ; run $NCPU

