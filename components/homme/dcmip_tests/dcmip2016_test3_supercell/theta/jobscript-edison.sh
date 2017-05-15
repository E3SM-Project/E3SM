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

# 4dg resolution
cp -f namelist-r400.nl input.nl
srun -n $NCPU $EXEC < input.nl
mv -f movies/dcmip2016_test31.nc movies/dcmip2016_test3_r400.nc

# 2dg resolution
cp -f namelist-r200.nl input.nl
srun -n $NCPU $EXEC < input.nl
mv -f movies/dcmip2016_test31.nc movies/dcmip2016_test3_r200.nc

# 1dg resolution
cp -f namelist-r100.nl input.nl
srun -n $NCPU $EXEC < input.nl
mv -f movies/dcmip2016_test31.nc movies/dcmip2016_test3_r100.nc

date
