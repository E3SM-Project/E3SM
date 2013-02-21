#!/bin/tcsh -f
#PBS -l nodes=16:ppn=8
#PBS -l walltime=1:00:00
#PBS -N swtc1
#PBS -j oe
#PBS -A FY081407
#PBX -a 1758
#XXX -W depend=afterany:jobid

#
#  Shallow water test case 6 "referece" test described in
#  homme/README 
#
set wdir = ~/scratch1/swtc6
set src = ~/codes/homme/build/sweqx
set input = ~/codes/homme/test/sw_conservative

set NCPU = 2
if ( ${?PBS_NODEFILE} ) then
   set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
endif
if ( ${?SLURM_NNODES} ) then
   # SLURM_NNODES  = number of nodes
   # hard to tell how many cores per nodes
   # set NCPU to zero, and mpirun will use the max allowed
   set NCPU = 0
endif
echo NCPU = $NCPU


set test_case = swtc6
set sub_case = 1 # default

#configure the model
cd $src
set configure = 1
if ( $configure ) then
  if (${?PNETCDF_PATH} ) then
     ./configure --with-netcdf=$NETCDF_PATH --with-pnetcdf=$PNETCDF_PATH NP=4 PLEV=1
  else
     ./configure --with-netcdf=$NETCDF_PATH NP=4 PLEV=1
  endif
  if ($status ) exit

  make clean
  make -j4 depends
  if ($status ) exit
endif

make -j4 sweqx
if ($status ) exit




mkdir $wdir
cd $wdir
mkdir movies

# defaults:
set smooth=0
set rk_stage=0
set nu = 0
set nu_s = -1  # defaults to nu
set LFTfreq = 0



set NE = 30
set nu = 1.0e15   
set hypervis_subcycle =  1

#leapfrog
# tstep=90  nu=0      NAN
#           nu=1e15   stable, but occasional 2dx noise at cube corners
#           nu=1.5e15 good
# tstep=80  nu=0      stable, but noisy
#           nu=1e15   good
#
set integration = explicit
set smooth = .05
set tstep = 80    
                  

#leapfrog-trap
#set integration = explicit
#set smooth = 0
#set tstep = 120
#set LFTfreq = 1

# RK2-m stage used by 3d code
#set integration = explicit
#set smooth = 0
#set tstep = 80
#set LFTfreq = 4

# RK-SSP
#set integration = runge_kutta
#set rk_stage = 3   
#set tstep = 30


set limiter = 0
set filter_freq = 0
set name = ${test_case}-NE${NE}-t${tstep}-limiter$limiter


set sfreq = 6
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`


sed s/ne=.\*/"ne = $NE"/  $input/swtc6high.nl |\
sed s/tstep.\*/"tstep = $tstep"/  |\
sed s/limiter_option.\*/"limiter_option = $limiter"/  |\
sed s/smooth.\*/"smooth = $smooth"/  |\
sed s/test_case.\*/"test_case = \'$test_case\' sub_case=$sub_case"/  |\
sed s/integration.\*/"integration = '$integration'"/  |\
sed s/rk_stage_user.\*/"rk_stage_user = $rk_stage"/  |\
sed s/LFTfreq.\*/"LFTfreq = $LFTfreq"/  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu_s"/  |\
sed s/filter_freq.\*/"filter_freq = $filter_freq"/  |\
sed s/hypervis_subcycle.\*/"hypervis_subcycle = $hypervis_subcycle"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> input.nl

date
mpirun -np $NCPU $src/sweqx < input.nl tee  sweq.out
date

