#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 3 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J d16-3-theta-r50    # job name
#SBATCH -o out_dcmip16-3.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 1152               # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

EXEC=../../../test_execs/theta-nlev40/theta-nlev40
NCPU=1152

date

hydrostatic="false"
#hydrostatic="true"

# 0.5 dg resolution
#sed -e "s/theta_hydrostatic_mode.*/theta_hydrostatic_mode=.${hydrostatic}./g" namelist-r50.nl >& input.nl
cp namelist-animation.nl input.nl

srun -n $NCPU $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r50.nc
mv HommeTime                  HommeTime_r50
mv max_w.pdf                  max_w_r50.pdf
mv max_precip.pdf             max_precip_r50.pdf
mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_r50.pdf
mv 5km_xsec.pdf               5km_xsec_r50.pdf
mv measurement_wmax.txt       measurement_wmax_r50.txt
mv measurement_time.txt       measurement_time_r50.txt
mv measurement_prect_rate.txt measurement_prect_rate_r50.txt

date

convert -density 100 -strip animation_5km_xsec.pdf tempimage-%03d.png
convert -layers OptimizePlus -delay 6 tempimage-*.png -loop 1 anim.gif
rm tempimage-*.png
