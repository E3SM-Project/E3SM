#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 2 on the NERSC Cori machine
#
# usage: sbatch jobscript-...

#SBATCH -J d16-2-theta        # job name
#SBATCH -o out_dcmip16-2.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 1280               # total number of mpi tasks requested
#SBATCH -p regular            # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00           # run time (hh:mm:ss)
#SBATCH -A m2618 # acme       # charge hours to this account
#SBATCH -C haswell            # use Haswell nodes
#SBATCH --qos=premium

EXEC=../../../test_execs/theta-nlev30/theta-nlev30
NCPU=1280

date

hydrostatic="false"
#hydrostatic="true"

# 1/2 dg resolution
sed -e "s/theta_hydrostatic_mode.*/theta_hydrostatic_mode=.${hydrostatic}./g" namelist-r50.nl >& input.nl
srun -n $NCPU $EXEC < input.nl

mv movies/dcmip2016_test21.nc movies/dcmip2016_test2_r50.nc
mv HommeTime HommeTime_r50

date
