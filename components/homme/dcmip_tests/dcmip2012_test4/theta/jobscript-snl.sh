#!/bin/tcsh 
#
#SBATCH --job-name dcmip4
#SBATCH --account=FY150001
#SBATCH -N 12
#SBATCH --time=0:10:00
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
#set namelist = nh-x1.nl           # 2592000 timesteps, 54 nodes - 7.4h
#set namelist = nh-x10.nl          #  648000 timsteps,  54 nodes - 115min
#set namelist = nh-x100.nl         #  64800  timesteps, 12 nodes - 51min 
set namelist = nh-x1000.nl        #  6480 timesteps, 12 nodes - 5.9min
#set namelist = h-x1.nl            #  4320 timesteps, 12 nodes, 3.3min

#set EXEC = ../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
#set namelist = preqx-x1.nl   # 4320 timesteps,   12 nodes, 2.2min


cp -f $namelist input.nl
date
mpirun -np $NCPU $EXEC < input.nl
date


ncl ps.ncl
ncl zeta.ncl

