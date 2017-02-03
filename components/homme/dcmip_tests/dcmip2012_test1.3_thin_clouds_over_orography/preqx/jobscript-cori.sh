#!/bin/bash
#
#   Jobscript for launching dcmip2012 test1-3 on the NERSC Cori machine
#
# usage: sbatch jobscript-...

#SBATCH -J dcmip1-3           # job name
#SBATCH -o out_dcmip1-3.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 768                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:05:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1
#SBATCH -C haswell            # use Haswell nodes

EXEC=../../../test_execs/preqx-nlev60-interp/preqx-nlev60-interp        # set name of executable
srun -n 768 $EXEC < ./namelist-default.nl                               # launch simulation

