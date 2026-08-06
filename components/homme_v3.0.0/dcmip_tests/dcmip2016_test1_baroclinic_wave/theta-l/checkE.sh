#!/bin/bash
#
#XXSBATCH --account=FY150001
#XXSBATCH -p ec
#SBATCH --account=condo
#SBATCH -p acme-centos6
#SBATCH -N 6
#SBATCH --time=2:00:00
#
# Anvil: 6 nodes, 1h for all runs
#
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=16M     #  Cori has 96GB per node. had to lower to 8M on 3K nodes
export MV2_ENABLE_AFFINITY=0
set | grep NODE


# hydrostatic
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30
case=nh
\rm -f input.nl
mkdir restart

function run { 
local tstep=$1

logfile1=log1-${tstep}.out
logfile2=log2-${tstep}.out
\rm -f $logfile1 $logfile2 

# run 15 days, write a restart file
namelist=checkE-run1.nl
sed s/tstep.\*1200/"tstep=$tstep"/  $namelist > input.nl
echo "tstep=" $tstep "run to 15 days and write restart file..."
\rm -f restart/*
mpirun $EXEC < input.nl  > $logfile1

grep SHA $logfile1

# restart from day 15, run 1h
let nmax=3*1200/$tstep
let nmax=2*3*1200/$tstep
rname=`ls -t restart/R*[0-9] | head -1 `


namelist=checkE-run2.nl
sed s/tstep.\*1200/"tstep=$tstep"/  $namelist |
sed s:nmax.\*:"nmax=$nmax":  |
sed s:restartfile.\*:"restartfile='$rname'":  > input.nl
echo "nmax=" $nmax " restartfile=" $rname " restart run..."
mpirun  $EXEC < input.nl > $logfile2

echo "diagnostics 3rd timsteps:"
grep "d/dt" $logfile2 | head -12 | tail -4
echo "diagnostics last timesteps:"
grep "d/dt" $logfile2 | tail -4
grep "E-E0" $logfile2 | tail -1

}

date
run 1200
run 600
run 300
run 150
run 75
run 36
run 18
date

