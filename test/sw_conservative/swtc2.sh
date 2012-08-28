#!/bin/tcsh -f
#XPBS -l nodes=100:ppn=4
#PBS -l nodes=200:ppn=2
#PBS -l walltime=8:00:00
#PBS -N swtc2
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


set test_case = swtc2
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
set nu_s = 0
set LFTfreq = 0

#ne=4 dt=600  

#set NE = 48
#set nu = 2e14
set NE = 120
set nu = 0
#set nu = 1e15      # needs dt < 150s
#set nu = 1.5e15   
#set nu = 2e15      # needs dt < 75s
#set nu = 3e15      # needs dt < 50s
#set nu = 5e15      # needs dt < 25s
set hypervis_subcycle =  1


set integration = explicit
#set smooth = .05
set smooth = 0 ; set LFTfreq = 1
#set tstep = 600
#set tstep = 300    # stable with nu=1e15
#set tstep = 180    # stable with nu=1e15
#set tstep = 90    # stable with nu=1e15
#set tstep = 45
set tstep = 20
#set tstep = 125
#set tstep = 70    # stable with nu=0

#set integration = runge_kutta
#set rk_stage = 2
#set tstep = 75

#set rk_stage = 3   
#set tstep = 30



set filter_freq = 0
set name = ${test_case}-NE${NE}-t${tstep}


set sfreq = 6
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`


sed s/ne=.\*/"ne = $NE"/  $input/swtc2.nl |\
sed s/tstep.\*/"tstep = $tstep"/  |\
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
> swtc2.nl

date
mpirun -np $NCPU $src/sweqx < swtc2.nl | tee  sweq.out
date

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out

mv -f swtc2.l1.errors  $name.l1.errors
mv -f swtc2.l2.errors  $name.l2.errors
mv -f swtc2.linf.errors  $name.linf.errors

tail -1 $name.l*.errors


