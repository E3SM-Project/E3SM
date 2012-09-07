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
#
set wdir = ~/scratch1/swtc6
set src = ~/codes/homme/build/sweqx
set input = ~/codes/homme/test/sw_conservative

set NCPU = 2
if ( ${?PBS_NODEFILE} ) then
   set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
endif
echo NCPU = $NCPU


set test_case = swtc6
set params = $input/Params4-1.inc    # 1 level


diff  $params  $src/../Params.inc
if ($status != 0) then
   echo "replacing Params.inc"
   cp $params $src/../Params.inc
endif



#configure the model
set configure = 1
if ( $configure ) then
  cd $src
  ./configure --with-netcdf=$NETCDF_PATH --with-pnetcdf=$PNETCDF_PATH NP=4 PLEV=1
  if ($status ) exit

  gmake clean
  gmake -j4 depends
  if ($status ) exit
endif

gmake -j2 sweqx
if ($status ) exit




mkdir $wdir
cd $wdir
mkdir movies

# defaults:
set smooth=0
set rk_stage=0
set nu = 0
set nu_s = 0
set LFTfreq = 0



set NE = 30
set nu = 1.5e15   
set hypervis_subcycle =  1


set integration = explicit
set smooth = .05
set tstep = 90    # stable with nu=1e15

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
sed s/test_case.\*/"test_case = \'$test_case\'"/  |\
sed s/integration.\*/"integration = '$integration'"/  |\
sed s/rk_stage_user.\*/"rk_stage_user = $rk_stage"/  |\
sed s/LFTfreq.\*/"LFTfreq = $LFTfreq"/  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu_s"/  |\
sed s/filter_freq.\*/"filter_freq = $filter_freq"/  |\
sed s/hypervis_subcycle.\*/"hypervis_subcycle = $hypervis_subcycle"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> swtc6.nl

date
mpirun -np $NCPU $src/sweqx < swtc6.nl | tee  sweq.out
date

