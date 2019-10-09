#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-0 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d20-preqx          # job name
#SBATCH -o out_dcmip2-0.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 64                 # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
srun -n 64 $EXEC < ./namelist-default.nl                               # launch simulation

