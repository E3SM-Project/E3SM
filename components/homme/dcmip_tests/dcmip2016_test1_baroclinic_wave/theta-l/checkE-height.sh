#!/bin/bash
#
#XXSBATCH --account=FY150001
#XXSBATCH -p ec
#SBATCH --account=condo
#SBATCH -p acme-centos6
#SBATCH -N 22
#SBATCH --time=2:00:00
#
# Anvil: 6 nodes, 1h for all runs
#
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=16M     #  Cori has 96GB per node. had to lower to 8M on 3K nodes
export MV2_ENABLE_AFFINITY=0
set | grep NODE


# 
#   run1: Run 15 days and write a restart file
#   run2: restart and run 2h in adiabatic mode for E convergence
#
#   this script uses dt=.5 and the explicit timestep, and was
#   used for the height coordinate tests (hcoord=1 in control_mod.F90)
#
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30
case=nh
\rm -f input.nl
mkdir restart

let  PPN=18
let  NCPU=384


if (( 1 )) ; then
  # run 15 days, write a restart file
  tstep=0.5
  logfile1=log1-${tstep}.out
  \rm -f $logfile1 

  namelist=checkE-run1-exp.nl
  sed s/tstep.\*1200/"tstep=$tstep"/  $namelist > input.nl
  echo "tstep=" $tstep "run to 15 days and write restart file..."
  \rm -f restart/*
  mpirun -bind-to=core -ppn $PPN -np $NCPU  $EXEC < input.nl  > $logfile1
  grep SHA $logfile1
  exit
fi


function run { 
local nmax=$1
local tstep=$2

logfile2=log2-$tstep.out
\rm -f  $logfile2 

rname=`ls -t restart/R*[0-9] | head -1 `

sfreq=1
if [ $nmax -gt 100 ]; then
sfreq=1440
fi


namelist=checkE-run2-exp.nl
sed s/tstep.\*0.5/"tstep=$tstep"/  $namelist |
sed s:nmax.\*:"nmax=$nmax":  |
sed s:statefreq.\*:"statefreq=$sfreq":  |
sed s:restartfile.\*:"restartfile='$rname'":  > input.nl
echo "nmax=" $nmax " restartfile=" $rname " restart run..."

mpirun -bind-to=core -ppn $PPN -np $NCPU    $EXEC < input.nl > $logfile2


echo "diagnostics 3rd timsteps: $logfile2"
grep "d/dt" $logfile2 | head -12 | tail -4
echo "diagnostics last timesteps:"
grep "d/dt" $logfile2 | tail -4
grep "E-E0" $logfile2 | tail -1

}

# .5,.25,.1,.05,.025, .01
date
run 14400  0.5
run 28800  0.25
run 57600  0.125
run 144000  0.050


#run 7200    0.5
#run 14400   0.25
#run 28800  0.125
#run 57600   0.050


date
