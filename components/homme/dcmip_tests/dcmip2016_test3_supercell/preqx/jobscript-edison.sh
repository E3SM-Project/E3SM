#!/bin/bash
#
#   Jobscript for launching dcmip2016 test3 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d16-3-preqx        # job name
#SBATCH -o out_dcmip16-3.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 192                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1
date

EXEC=../../../test_execs/preqx-nlev40-interp/preqx-nlev40-interp        # set name of executable
cp ./namelist-default.nl input.nl
srun -n 192 $EXEC < input.nl                                            # launch simulation

date
