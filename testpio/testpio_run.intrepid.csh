#!/bin/csh -f

# --- user can set these, the rest are derived ---
set ntasks = 16
set mthrds = 1
set project = CCESDev
set input_file = "testpio_in.pb02"
# ----------------------

set max_tasks_per_node = 4
@ nodes = (${ntasks} * ${mthrds}) / ${max_tasks_per_node}
if ((${nodes} * ${max_tasks_per_node}) < (${ntasks} * ${mthrds})) then
  @ nodes = ${nodes} + 1
endif

if (${mthrds} == 1) then
  set mode = vn
else if (${mthrds} == 2) then
  set mode = dual
else if (${mthrds} == 4) then
  set mode = smp
else
  echo "ERROR illegal max thread count ${mthrds}"
  exit 1
endif

echo "case = $0"
set jobsub = `qstat --header User:Command | grep $0 | wc -l`
set LID = "`date +%y%m%d-%H%M%S`"
if (${jobsub} == 0) then
    echo "qsub -n ${nodes} -t 30 -q prod-devel -A $project --mode script $0"
    qsub -n ${nodes} -t 30 -q prod-devel -A $project --mode script -o testpio.out.$LID $0
    exit 0
endif

# ----------------------

set srcdir = `pwd`
set wrkdir = "~/scratch/testpio"

if (! -d $wrkdir) mkdir -p $wrkdir
cd $wrkdir
rm -f ./testpio
cp -f $srcdir/testpio ./testpio
rm -f ./testpio_in
cp -f $srcdir/${input_file} ./testpio_in
if (! -d none) mkdir none
rm -r -f none/*

cobalt-mpirun -np ${ntasks} -mode ${mode} -verbose 2 -cwd `pwd` -env "XLSMPOPTS=stack=64000000 OMP_NUM_THREADS=${mthrds} DCMF_COLLECTIVES=1 BG_MAPPING=TXYZ" ./testpio

cp testpio.out.$LID $srcdir/

