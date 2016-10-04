#!/bin/bash
#
#   Jobscript for launching dcmip2012 test3-1 on the NERSC Edison machine
#
#SBATCH -J dcmip3-1           # job name
#SBATCH -o out_dcmip3-1.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 1344               # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

EXEC=../../test_execs/preqx-nlev20/preqx-nlev20                         # set name of executable
srun -n 1344 $EXEC < ./namelist-default.nl                              # launch simulation

