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

set test_case = swtc2

#configure the model
cd $src
set configure = 0
if ( $configure ) then
  cd $src
  ./configure --with-netcdf=$NETCDF_PATH --with-pnetcdf=$PNETCDF_PATH NP=4 PLEV=1
  if ($status ) exit

  make clean
  make -j4 depends
  if ($status ) exit
endif

make -j2 sweqx
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

set hypervis_subcycle =  1
set integration = explicit

set NE = 20
#set nu = 5.1e15
set tstep = 144

### leapfrog
#set smooth = 0.05 ; set LFTfreq = 0

### leapfrog-trapazoidal 
set smooth = 0 ; set LFTfreq = 1

### RK (note: needs viscosity to be stable)
#set integration = runge_kutta
#set rk_stage = 2
#set tstep = 75

set filter_freq = 0
set name = ${test_case}-NE${NE}-t${tstep}


set sfreq = 6
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`


sed s/ne=.\*/"ne = $NE"/  $input/swtc2.nl |\
sed s/tstep.\*/"tstep = $tstep"/  |\
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
mpirun -np $NCPU $src/sweqx < input.nl | tee  sweq.out
date

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out

mv -f swtc2.l1.errors  $name.l1.errors
mv -f swtc2.l2.errors  $name.l2.errors
mv -f swtc2.linf.errors  $name.linf.errors

tail -1 $name.l*.errors


