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
set src = ~/codes/homme/build/sweqx
set input = ~/codes/homme/test/sw_conservative

set NCPU = 0
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

if ( $build == 1 ) then
   cd $src
   ./configure --enable-blas --enable-lapack --with-netcdf=$NETCDF_PATH \
    --with-pnetcdf=$PNETCDF_PATH NP=4 PLEV=1   --enable-energy-diagnostics
   make depends
   make clean
   make -j4 sweqx
   exit
endif
if ( $make == 1 ) then
   cd src
   make -j4 sweqx
   if ($status) exit
endif


mkdir $wdir
cd $wdir
mkdir movies

# defaults:
set smooth=0
set LFTfreq = 0
set NE = 30
set ndays = 2 # official test case length 15 days
set nu = 1.5e15   
set nu_s = $nu
set NPDG = 0

# turn on hybrid CG/DG  
set NPDG = 3  # 0 = disabled
set nu_s = 0  # dont use viscsosity with DG code
set nu = 10e15   



set hypervis_subcycle =  1
set integration = runge_kutta
set rk_stage = 3
set tstep = 30

set limiter = 0
set filter_freq = 0
set name = ${test_case}-NE${NE}-t${tstep}


set sfreq = 1
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`

# output units: 0,1,2 = timesteps, days, hours
set OUTUNITS = 2  
set OUTFREQ =  4

sed s/ne=.\*/"ne = $NE  npdg=$NPDG"/  $input/swtc6high.nl |\
sed s/tstep.\*/"tstep = $tstep"/  |\
sed s/limiter_option.\*/"limiter_option = $limiter"/  |\
sed s/smooth.\*/"smooth = $smooth"/  |\
sed s/test_case.\*/"test_case = \'$test_case\'"/  |\
sed s/ndays.\*/"ndays = $ndays"/  |\
sed s/integration.\*/"integration = '$integration'"/  |\
sed s/rk_stage_user.\*/"rk_stage_user = $rk_stage"/  |\
sed s/LFTfreq.\*/"LFTfreq = $LFTfreq"/  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu_s"/  |\
sed s/filter_freq.\*/"filter_freq = $filter_freq"/  |\
sed s/hypervis_subcycle.\*/"hypervis_subcycle = $hypervis_subcycle"/  |\
sed s/output_frequency.\*/"output_frequency = $OUTFREQ"/  |\
sed s/output_timeunits.\*/"output_timeunits = $OUTUNITS"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> swtc5.nl

date
rm -f movies/swtc5?.nl
mpirun -np $NCPU $src/sweqx < swtc5.nl | tee  sweq.out
date


