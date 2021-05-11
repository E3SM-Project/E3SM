#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 2 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d16-2-preqx        # job name
#SBATCH -o out_dcmip16-2.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 640                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge account 

EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
NCPU=640

date

# 1dg resolution
cp namelist-r100.nl input.nl
srun -n $NCPU $EXEC < input.nl

mv movies/dcmip2016_test21.nc movies/dcmip2016_test2_r100.nc
mv HommeTime HommeTime_r100

date
