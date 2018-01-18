#!/bin/tcsh 
#
#SBATCH --job-name d16-3-theta
#SBATCH --account=FY150001
#SBATCH -N 36
#SBATCH --time=0:30:00
#SBATCH -p ec
#PBS -q acme
#PBS -l walltime=60:00                                                                                           
#PBS -l nodes=25    

set OMP_NUM_THREADS = 1
set NCPU = 640
if ( ${?PBS_ENVIRONMENT} ) then   # anvil
  set NCPU = $PBS_NNODES
  if ( $PBS_ENVIRONMENT == PBS_BATCH ) cd $PBS_O_WORKDIR     
endif
if ( ${?SLURM_NNODES} ) then   # redsky
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

set EXEC = ../../../test_execs/theta-nlev40/theta-nlev40

date

# 4dg resolution
\cp -f namelist-r400.nl input.nl
mpirun -np 320 $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

\mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r400.nc
\mv HommeTime                  HommeTime_r400
\mv max_w.pdf                  max_w_r400.pdf
\mv max_precip.pdf             max_precip_r400.pdf
\mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_r400.pdf
\mv 5km_xsec.pdf               5km_xsec_r400.pdf
\mv measurement_wmax.txt       measurement_wmax_r400.txt
\mv measurement_time.txt       measurement_time_r400.txt
\mv measurement_prect_rate.txt measurement_prect_rate_r400.txt

date

# 2dg resolution
\cp -f namelist-r200.nl input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

\mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r200.nc
\mv HommeTime                  HommeTime_r200
\mv max_w.pdf                  max_w_r200.pdf
\mv max_precip.pdf             max_precip_r200.pdf
\mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_r200.pdf
\mv 5km_xsec.pdf               5km_xsec_r200.pdf
\mv measurement_wmax.txt       measurement_wmax_r200.txt
\mv measurement_time.txt       measurement_time_r200.txt
\mv measurement_prect_rate.txt measurement_prect_rate_r200.txt

date

# 1dg resolution
\cp -f namelist-r100.nl input.nl
mpirun -np $NCPU $EXEC < input.nl
ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

\mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_r100.nc
\mv HommeTime                  HommeTime_r100
\mv max_w.pdf                  max_w_r100.pdf
\mv max_precip.pdf             max_precip_r100.pdf
\mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_r100.pdf
\mv 5km_xsec.pdf               5km_xsec_r100.pdf
\mv measurement_wmax.txt       measurement_wmax_r100.txt
\mv measurement_time.txt       measurement_time_r100.txt
\mv measurement_prect_rate.txt measurement_prect_rate_r100.txt

date

