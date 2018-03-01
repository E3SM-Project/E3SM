#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-2 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d22-theta          # job name
#SBATCH -o out_dcmip2-2.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 320                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:20:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

EXEC=../../../test_execs/theta-nlev60/theta-nlev60
NCPU=320

date

# hydrostatic theta
cp  -f namelist-h.nl input.nl
srun -n $NCPU $EXEC < input.nl
ncl  plot_lon_vs_z.ncl
mv  -f movies/dcmip2012_test2_21.nc  movies/hydro_dcmip2012_test2_21.nc
mv  -f dcmip2012_test2_2_T_t10.pdf   hydro_T_t10.pdf

# nonhydrostatic theta
cp -f namelist-nh.nl input.nl
srun -n $NCPU $EXEC < input.nl
ncl plot_lon_vs_z.ncl
mv -f movies/dcmip2012_test2_21.nc  movies/nonhydro_dcmip2012_test2_21.nc
mv -f dcmip2012_test2_2_T_t10.pdf   nonhydro_T_t10.pdf

date
