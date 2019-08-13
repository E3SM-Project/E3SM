#!/bin/bash
#
#SBATCH --job-name dcmip2016-3
#SBATCH -N 512
#SBATCH -C knl
#SBATCH --time=4:00:00
#SBATCH -q regular
#
#                   init time, run time, I/O per snapshot
#  NE=900  namelist-r100-x4.nl   512 nodes:  530s,  1h,  450s                       2h-run: 1.3h
#  NE=1800 namelist-r100-x2.nl   512 nodes:  0.4h, 1.9h, 891s  (geo,w,Q2,Q3 only)   2h-run: 3.5h
#  NE=3600 namelist-r100-x1.nl   1536 nodes: 1.1h, 5.1h, 1.1h    (geo,w,Q2,Q3 only)    12h?
#
#

#
#  mpi run command
#
export OMP_STACKSIZE=16M     #  Cori has 96GB per node. had to lower to 8M on 3K nodes
export OMP_NUM_THREADS=2
PER_NODE=64          #  MPI per node

# number of virtual cores per MPI task
VC_PER_MPI=256  # Set this to 272 if using PER_NODE divides 272 instead of 256
let VC_PER_MPI/=$PER_NODE

export KMP_AFFINITY="granularity=core,scatter"
bind=--cpu_bind=core

# compute number of MPI tasks
if [ -n "$SLURM_NNODES" ]; then
  NNODES=$SLURM_NNODES
else
  NNODES=1
fi


NMPI=$NNODES
let NMPI*=$PER_NODE

# default DCMIP executable uses interpolated output
# need to compile different version with native grid output
# assume wdir already configured for Cori, but reconfigure theta-l
# target with correct levels and native grid output:
if (( 0 )) ; then
  echo running cmake
  cd ../../..
  cmake -DQSIZE_D=3 -DPREQX_PLEV=40 -DPREQX_NP=4  \
   -DBUILD_HOMME_SWEQX=FALSE  -DPREQX_USE_PIO=TRUE             \
   -DPREQX_USE_ENERGY=FALSE  .
  make -j4 theta-l
  cd dcmip_tests/dcmip2016_test3_supercell/theta-l
  exit
fi


EXEC=../../../src/theta-l/theta-l


mpirun="srun -n $NMPI -N $NNODES -c $VC_PER_MPI $bind"
echo mpi commnand:
echo $mpirun

namelist=namelist-r100-x4.nl
\cp -f $namelist input.nl
date
$mpirun  $EXEC < input.nl
date

