#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 1 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d16-3-theta          # job name
#SBATCH -o out_dcmip16-3.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 320                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:20:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

EXEC=../../../test_execs/theta-nlev40/theta-nlev40
NCPU=320

date

# hydrostatic theta
cp  -f namelist-h.nl input.nl
srun -n $NCPU $EXEC < input.nl
#ncl  plot_lon_vs_z.ncl
#mv  -f movies/dcmip2012_test2_11.nc  movies/hydro_dcmip2012_test2_11.nc
#mv  -f dcmip2012_test2_1_T_t10.pdf   hydro_T_t10.pdf

# nonhydrostatic theta
cp -f namelist-nh.nl input.nl
srun -n $NCPU $EXEC < input.nl
#ncl plot_lon_vs_z.ncl
#mv -f movies/dcmip2012_test2_11.nc  movies/nonhydro_dcmip2012_test2_11.nc
#mv -f dcmip2012_test2_1_T_t10.pdf   nonhydro_T_t10.pdf

date
