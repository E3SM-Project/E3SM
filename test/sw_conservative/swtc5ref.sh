#!/bin/tcsh -f
#XPBS -l nodes=100:ppn=4
#PBS -l nodes=200:ppn=2
#PBS -l walltime=1:00:00
#PBS -N swtc1
#PBS -j oe
#PBS -A FY081407
#PBX -a 1758
#XXX -W depend=afterany:jobid

#
#  Shallow water test case 6 "referece" test described in
#  homme/README 
#  Mark Taylor 2010/10
#
set wdir = ~/scratch1/swtc5
set HOMME = ~/codes/homme
set MACH = $HOMME/cmake/machineFiles/redsky.cmake
set src = $HOMME/build/sweqx
set input = $HOMME/test/sw_conservative
mkdir $wdir
cd $wdir
mkdir movies

set NCPU = 32
if ( ${?PBS_NODEFILE} ) then
   set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
endif
echo NCPU = $NCPU


set test_case = swtc5


set build = 0
set make = 1
if ( $#argv >= 1) then
  if ( $1 == 'build' ) set build = 1
endif

# configure build system:
# if ( $build == 1 ) then
#    cd $src
#    ./configure --enable-blas --enable-lapack --with-netcdf=$NETCDF_PATH \
#      --with-pnetcdf=$PNETCDF_PATH NP=4 PLEV=1   --enable-energy-diagnostics
#    make depends
#    make clean
#    make -j4 sweqx
#    exit
# endif
# if ( $make == 1 ) then
#    cd src
#    make -j4 sweqx
#    if ($status) exit
# endif
# set exe = $src/sweqx
#cmake:
cd $wdir
if ( $build == 1 ) then
   rm -rf CMakeFiles CMakeCache.txt
   cmake -C $MACH -DSWEQX_PLEV=1  -DSWEQX_NP=4 $HOMME
   exit
endif
if ( $make == 1 ) then
   make -j4 sweqx
    if ($status) exit
endif
set exe = $wdir/src/sweqx/sweqx




# defaults:
set smooth=0
set rk_stage=0
set nu = 0
set nu_s = -1  # defaults to nu
set LFTfreq = 0
set hypervis_subcycle =  1

set NE = 30

#LF
set tstep = 90
set nu = 1.5e15   
set integration = explicit
set smooth = .05

# rk
#set smooth=0
#set integration = runge_kutta
#set rk_stage = 2
#set tstep = 75
#set rk_stage = 3   
#set tstep = 30
#set nu=0



set limiter = 0
set filter_freq = 0
set name = ${test_case}-NE${NE}-t${tstep}


set sfreq = 6
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`


sed s/ne=.\*/"ne = $NE"/  $input/swtc6high.nl |\
sed s/tstep.\*/"tstep = $tstep"/  |\
sed s/limiter_option.\*/"limiter_option = $limiter"/  |\
sed s/smooth.\*/"smooth = $smooth"/  |\
sed s/test_case.\*/"test_case = \'$test_case\'"/  |\
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
rm -f movies/swtc5?.nc
mpirun -np $NCPU $exe < input.nl | tee  sweq.out
date


