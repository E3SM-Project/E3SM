#!/bin/tcsh -f
#XPBS -l nodes=100:ppn=4
#PBS -l nodes=200:ppn=2
#PBS -l walltime=8:00:00
#PBS -N swtc1
#PBS -j oe
#PBS -A FY081407
#PBX -a 1758
#XXX -W depend=afterany:jobid

#
#  Shallow water test case 6 Runge-Kutta test
#  used to check stability of m-stage RK schemes
#
set wdir = ~/scratch1/sweqx
set src = ~/codes/homme/build.Linux
set input = ~/codes/homme/test/sw_conservative
set NCPU = 2

if ( ${?SNLSYSTEM} ) then
if ( $SNLSYSTEM == rose | $SNLSYSTEM == thunderbird ) then
   if ( ${?PBS_NODEFILE} ) then
      set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
   endif
endif
endif
echo NCPU = $NCPU


set test_case = swtc6
set params = $input/Params4-1.inc    # 1 level


diff  $params  $src/../Params.inc
if ($status != 0) then
   echo "replacing Params.inc"
   cp $params $src/../Params.inc
endif



cd $src
rm -f sweqx
make -j2 sweqx
mkdir $wdir
cd $wdir
mkdir movies

# defaults:
set smooth=0
set rk_stage=0
set nu = 0
set nu_s = -1     # defaults to nu
set LFTfreq = 0



#set NE = 48
#set nu = 2e14
set NE = 30
#set nu = 0
set nu = 1e15      # needs dt < 150s
#set nu = 1.5e15   
#set nu = 2e15      # needs dt < 75s
#set nu = 3e15      # needs dt < 50s
#set nu = 5e15      # needs dt < 25s
set hypervis_subcycle =  1


set integration = explicit
#set smooth = .05
set smooth = 0 ; set LFTfreq = 1
set tstep = 90    # stable with nu=1e15
#set tstep = 125
#set tstep = 70    # stable with nu=0

#set integration = runge_kutta
#set rk_stage = 2
#set tstep = 75

#set rk_stage = 3   
#set tstep = 30



set limiter = 0

#limiter_option=1  adv_elem_monotone  hyper_zero
#limiter_option=2  adv_zero           hyper_elem_monotone
#limiter_option=3  adv_elem_monotone  hyper_elem_monotone
#limiter_option=4  adv_zero           hyper_zero
#limiter_option=5  adv_monotone       hyper_zero


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
sed s/integration.\*/"integration = $integration"/  |\
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

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out

