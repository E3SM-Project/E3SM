#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 3 on the NERSC Edison machine
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

hydrostatic="false"
#hydrostatic="true"

# 4dg resolution
sed -e "s/theta_hydrostatic_mode.*/theta_hydrostatic_mode=.${hydrostatic}./g" namelist-r400.nl >& input.nl
srun -n 320 $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r400.nc
mv HommeTime                  HommeTime_r400
mv max_w.pdf                  max_w_r400.pdf
mv max_precip.pdf             max_precip_r400.pdf
mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_r400.pdf
mv 5km_xsec.pdf               5km_xsec_r400.pdf
mv measurement_wmax.txt       measurement_wmax_r400.txt
mv measurement_time.txt       measurement_time_r400.txt
mv measurement_prect_rate.txt measurement_prect_rate_r400.txt

date

# 2dg resolution
sed -e "s/theta_hydrostatic_mode.*/theta_hydrostatic_mode=.${hydrostatic}./g" namelist-r200.nl >& input.nl
srun -n $NCPU $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r200.nc
mv HommeTime                  HommeTime_r200
mv max_w.pdf                  max_w_r200.pdf
mv max_precip.pdf             max_precip_r200.pdf
mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_r200.pdf
mv 5km_xsec.pdf               5km_xsec_r200.pdf
mv measurement_wmax.txt       measurement_wmax_r200.txt
mv measurement_time.txt       measurement_time_r200.txt
mv measurement_prect_rate.txt measurement_prect_rate_r200.txt

date

# 1dg resolution
sed -e "s/theta_hydrostatic_mode.*/theta_hydrostatic_mode=.${hydrostatic}./g" namelist-r100.nl >& input.nl
srun -n $NCPU $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r100.nc
mv HommeTime                  HommeTime_r100
mv max_w.pdf                  max_w_r100.pdf
mv max_precip.pdf             max_precip_r100.pdf
mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_r100.pdf
mv 5km_xsec.pdf               5km_xsec_r100.pdf
mv measurement_wmax.txt       measurement_wmax_r100.txt
mv measurement_time.txt       measurement_time_r100.txt
mv measurement_prect_rate.txt measurement_prect_rate_r100.txt

date
