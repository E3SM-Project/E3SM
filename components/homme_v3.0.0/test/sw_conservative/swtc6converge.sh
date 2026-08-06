#!/bin/tcsh -f
#PBS -l nodes=80:ppn=4
#XPBS -l nodes=200:ppn=2
#PBS -l walltime=4:00:00
#PBS -N swtc1
#PBS -j oe
#PBS -A FY081407
#PBX -a 1758
#XXX -W depend=afterany:jobid

#
#  Shallow water test case 6 script used for convergence study
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


#set test_case = swtc6
set test_case = swtc5
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
set nu_s = -1  # defaults to nu
set LFTfreq = 0
set hypervis_subcycle =  1

# T85 is like NE=30
# T170 is like NE=60
# T213 is like NE=75

# swtc5  E-E0/E0 NE=8     NE=15
#  640    0.7234E-06
#  320   -0.2386E-06      -0.2351E-06 
#  160   -0.4788E-07      -0.4651E-07 
#   80   -0.8451E-08      -0.8731E-08            
#   40   -0.1230E-08      -0.1567E-08  
#   20                    -0.2370E-09  3min            
#   10                    -0.3364E-10
#    5                    -0.5045E-11
#    2                    -0.4834E-12
#    1                    -0.9498E-13            
#
#  ENS -ENS0/ENS0  NE=15    note: initial value 7.2046e-13 
#  320   0.1063E-03
#  160   0.1063E-03
#   80   0.1063E-03
#   40   0.1063E-03
#
 
#set NE = 30
#set tstep = 80

set NE = 15
set tstep = 5

#set NE = 8
#set tstep = 40
#set NE = 4
#set tstep = 600

#set NE = 30
#set nu = 0
#set nu = 1e15      # needs dt < 150s
#set nu = 1.5e15   
#set nu = 2e15      # needs dt < 75s
#set nu = 3e15      # needs dt < 50s
#set nu = 5e15      # needs dt < 25s
#set hypervis_subcycle =  1
#set tstep = 90    # stable with nu=1e15
#set tstep = 125
#set tstep = 70    # stable with nu=0


#set NE = 48
#set nu = 2e14





#
#  time step choice
#
set integration = explicit
#set smooth = .05
set LFTfreq = 1; set smooth = .00


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
set name = ${test_case}-NE${NE}-t${tstep}
if ( $LFTfreq == 1 ) then
   set name = ${test_case}-NE${NE}-t${tstep}-LFT
endif

set sfreq = 8
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
sed s/nlat.\*/"nlat = 256"/  |\
sed s/nlon.\*/"nlon = 512"/  |\
sed s/output_end_time.\*/"output_end_time = 15"/  |\
sed s/filter_freq.\*/"filter_freq = $filter_freq"/  |\
sed s/hypervis_subcycle.\*/"hypervis_subcycle = $hypervis_subcycle"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> swtc6.nl

date
mpirun -np $NCPU $src/sweqx < swtc6.nl | tee  sweq.out
aprun -n $NCPU $src/sweqx < swtc6.nl | tee  sweq.out
date

mv -f movies/${test_case}1.nc  $name.nc

mv -f sweq.out $name.out
mv -f sweq.mass $name.mass
mv -f sweq.ke $name.ke
mv -f sweq.pe $name.pe
mv -f sweq.penst $name.penst
mv -f sweq.pv $name.pv
mv -f sweq.div $name.div




