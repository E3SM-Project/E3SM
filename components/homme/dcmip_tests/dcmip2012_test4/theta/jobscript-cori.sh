#!/bin/bash
#
#   Jobscript for launching dcmip2012 test3-1 on the NERSC Cori machine
#
# usage: sbatch jobscript-...

#SBATCH -J dcmip2-0           # job name
#SBATCH -o out_dcmip2-0.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 600                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:05:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1
#SBATCH -C haswell            # use Haswell nodes

#EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
srun -n 600 $EXEC < ./namelist-default.nl                               # launch simulation

