#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 1 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d16-3-theta        # job name
#SBATCH -o out_dcmip16-3.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 640                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

EXEC=../../../test_execs/theta-nlev40/theta-nlev40
NCPU=640

date

# hydrostatic theta
#cp  -f namelist-h.nl input.nl
#srun -n $NCPU $EXEC < input.nl

# nonhydrostatic theta
cp -f namelist-nh.nl input.nl
srun -n $NCPU $EXEC < input.nl

date
