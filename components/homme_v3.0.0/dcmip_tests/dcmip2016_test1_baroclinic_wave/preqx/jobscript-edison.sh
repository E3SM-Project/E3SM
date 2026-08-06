#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 1 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d16-1-preqx        # job name
#SBATCH -o out_dcmip16-1.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 640                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:20:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
NCPU=640

date

cp -f namelist-r100.nl input.nl
srun -n $NCPU $EXEC < input.nl

date
