#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-0 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d20-theta          # job name
#SBATCH -o out_dcmip2-0.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 432                # total number of mpi tasks requested
#SBATCH -p regular            # queue (partition)
#SBATCH -t 03:00:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

#SXXBATCH -n 48               # total number of mpi tasks requested
#SXXBATCH -t 00:05:00         # run time (hh:mm:ss)

# hydrostatic theta
date
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
cp ./namelist-h.nl input.nl
srun -n 432 $EXEC < ./input.nl
ncl plot_z_lon.ncl
ncl test200-range.ncl
mv -f dcmip2012_test2_0_u_t6.00.pdf     hydro_test2_0_u_z.pdf
mv -f movies/dcmip2012_test2_01.nc.pdf  hydro_test2_0_u.pdf
mv -f movies/dcmip2012_test2_01.nc      movies/hydro_dcmip2012_test2_01.nc

# nonhydrostatic theta
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
cp ./namelist-nh.nl input.nl
srun -n 432 $EXEC < ./input.nl
date
ncl plot_z_lon.ncl
ncl test200-range.ncl
mv -f dcmip2012_test2_0_u_t6.00.pdf     nonhydro_test2_0_u_z.pdf
mv -f movies/dcmip2012_test2_01.nc.pdf  nonhydro_test2_0_u.pdf
mv -f movies/dcmip2012_test2_01.nc      movies/nonhydro_dcmip2012_test2_01.nc
date

