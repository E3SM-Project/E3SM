#!/bin/bash
#
#   Jobscript for launching dcmip2012 test3-1 on the NERSC Cori machine
#
#SBATCH -J d31-theta          # job name
#SBATCH -o out_dcmip3-1.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 64                 # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1
#SBATCH -C haswell            # use Haswell nodes

EXEC=../../../test_execs/theta-nlev20/theta-nlev20
NCPU=64

#############################################################################
# theta (hydrostatic
#############################################################################
namelist=./namelist-h.nl
cp -f $namelist input.nl
srun -n $NCPU $EXEC < input.nl

ncl plot_omega.ncl
ncl plot_theta.ncl

mv -f test31_omega.pdf                 hydro_test31_omega.pdf
mv -f dcmip2012_test3_theta_diff.pdf   hydro_test3_theta_diff.pdf
mv -f movies/dcmip2012_test31.nc       movies/hydro_dcmip2012_test31.nc

#############################################################################
# theta-nh
#############################################################################
namelist=./namelist-nh.nl
cp -f $namelist input.nl
srun -n $NCPU $EXEC < input.nl

ncl plot_omega.ncl
ncl plot_theta.ncl

mv -f test31_omega.pdf                 nonhydro_test31_omega.pdf
mv -f dcmip2012_test3_theta_diff.pdf   nonhydro_test3_theta_diff.pdf
mv -f movies/dcmip2012_test31.nc        movies/nonhydro_dcmip2012_test31.nc

