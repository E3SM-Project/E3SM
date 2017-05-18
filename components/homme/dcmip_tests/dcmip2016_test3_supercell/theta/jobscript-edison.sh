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
#cp -f namelist-r400.nl input.nl
#srun -n 384 $EXEC < input.nl
#mv -f movies/dcmip2016_test31.nc movies/dcmip2016_test3_r400.nc
#ncl plot_supercell_wvel.ncl

#mv HommeTime HommeTime_r400
#mv measurement_wmax.pdf wmax_r400.pdf
#mv measurement_wmax.txt measurement_wmax_r400.txt
#mv measurement_time.txt measurement_time_r400.txt
#mv measurement_prect_rate.txt measurement_prect_rate.txt_r400.txt

# 2dg resolution
cp -f namelist-r200.nl input.nl
srun -n $NCPU $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_5km_xsec.ncl

mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r200.nc
mv HommeTime HommeTime_r200
mv measurement_wmax.pdf wmax_r200.pdf
mv 5km_xsec.pdf 5km_xsec_r200.pdf
mv measurement_wmax.txt measurement_wmax_r200.txt
mv measurement_time.txt measurement_time_r200.txt
mv measurement_prect_rate.txt measurement_prect_rate_r200.txt

# 1dg resolution
cp -f namelist-r100.nl input.nl
srun -n $NCPU $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_5km_xsec.ncl

mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r100.nc
mv HommeTime HommeTime_r100
mv measurement_wmax.pdf wmax_r100.pdf
mv 5km_xsec.pdf 5km_xsec_r100.pdf
mv measurement_wmax.txt measurement_wmax_r100.txt
mv measurement_time.txt measurement_time_r100.txt
mv measurement_prect_rate.txt measurement_prect_rate_r100.txtt


date
