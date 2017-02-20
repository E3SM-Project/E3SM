#!/bin/tcsh 
#
#SBATCH --job-name d4
#SBATCH --account=FY150001
#SBATCH -N 54
#SBATCH --time=00:05:00
#SBATCH -p ec


set OMP_NUM_THREADS = 1
set NCPU = 40 
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

# 2 min on 54 nodes on skybridge
# NH model x1000 case
set EXEC = ../../../test_execs/theta-nlev30/theta-nlev30
set namelist = nh-x1000.nl

# 320? minutes on 54 nodes on skybridge
# NH model x1 case
#set namelist = nh-x1.nl

# 80 min on 54 nodes on skybridge
# NH model x10 case
#set namelist = nh-x10.nl

# 15 min on 54 nodes on skybridge
# NH model x100 case
#set namelist = nh-x100.nl

# <10 minutes on 12 nodes on skybridge
# H  model x1 case
#set namelist = h-x1.nl


cp -f $namelist input.nl
date
mpirun -np $NCPU $EXEC < input.nl
date


cd movies
ncl ../ps.ncl
