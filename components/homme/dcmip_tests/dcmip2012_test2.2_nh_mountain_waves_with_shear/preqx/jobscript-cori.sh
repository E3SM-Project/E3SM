#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-2 on the NERSC Cori machine
#
# usage: sbatch jobscript-...

#SBATCH -J d22-preqx          # job name
#SBATCH -o out_dcmip2-2.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 192                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1
#SBATCH -C haswell            # use Haswell nodes
date

EXEC=../../../test_execs/preqx-nlev60-interp/preqx-nlev60-interp        # set name of executable
cp ./namelist-default.nl input.nl
srun -n 192 $EXEC < input.nl                                            # launch simulation
ncl plot_lon_vs_z.ncl

date
