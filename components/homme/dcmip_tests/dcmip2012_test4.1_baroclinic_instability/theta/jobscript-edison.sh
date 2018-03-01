#!/bin/bash
#
#   Jobscript for launching dcmip2012 test4-1 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d41-theta          # job name
#SBATCH -o out_dcmip4-1.o%j   # output and error file name (%j expands to jobID)
#SBATCH -n 192                # total number of mpi tasks requested
#SBATCH -p regular            # queue (partition) -- normal, development, etc.
#SBATCH -t 01:30:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

NCPU=192
date

#  preqx
EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
namelist=preqx-x1.nl   # 4320 timesteps,   12 nodes, 2.2min
cp -f $namelist input.nl
srun -n $NCPU $EXEC < input.nl
ncl ps.ncl   'fname="movies/preqx-X1-dcmip2012_test41.nc"'
ncl zeta.ncl 'fname="movies/preqx-X1-dcmip2012_test41.nc"'
mv -f zeta.pdf preqx-X1-zeta.pdf
mv -f ps.pdf preqx-X1-ps.pdf

date
# hydrostatic
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
namelist=h-x1.nl            #  4320 timesteps, 12 nodes, 3.3min?
cp -f $namelist input.nl
srun -n $NCPU $EXEC < input.nl
ncl ps.ncl   'fname="movies/hydro-X1-dcmip2012_test41.nc"'
ncl zeta.ncl 'fname="movies/hydro-X1-dcmip2012_test41.nc"'
mv -f zeta.pdf hydro-X1-zeta.pdf
mv -f ps.pdf hydro-X1-ps.pdf

date
# nonhydro X1000
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
namelist=nh-x1000.nl        #  6480 timesteps, 12 nodes - 5.9min?
cp -f $namelist input.nl
srun -n $NCPU $EXEC < input.nl
ncl ps.ncl   'fname="movies/nonhydro-X1000-dcmip2012_test41.nc"'
ncl zeta.ncl 'fname="movies/nonhydro-X1000-dcmip2012_test41.nc"'
mv -f zeta.pdf nonhydro-X1000-zeta.pdf
mv -f ps.pdf nonhydro-X1000-ps.pdf

date
# nonhydro X100
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
namelist=nh-x100.nl        #  64800  timesteps, 12 nodes - 51min?
cp -f $namelist input.nl
srun -n $NCPU $EXEC < input.nl
ncl ps.ncl   'fname="movies/nonhydro-X100-dcmip2012_test41.nc"'
ncl zeta.ncl 'fname="movies/nonhydro-X100-dcmip2012_test41.nc"'
mv -f zeta.pdf nonhydro-X100-zeta.pdf
mv -f ps.pdf nonhydro-X100-ps.pdf

date
