#!/bin/csh -f

source global_variables.csh

setenv case t${casestr}1_${NTHRDS_VAL}
setenv casebase ${mach}_${compset}_${res}
setenv case_dir ${casedir}/${case}_${casebase}

echo "Finding balanced layouts for these total task counts: "$TARGET_TASKS
echo "Results will be written to: "$results_dir

set numFiles = `ls -lt /glade/u/home/mickelso/test/restults/*.dat | wc -l`
if ($numFiles > 2) then
  set curDir = $PWD
  echo "Cleaning *.dat files from the results directory: "$results_dir
  cd $results_dir
  rm -fr *.dat 
  cd $curDir
endif

echo "Finding new layouts ..."
if ($DYCORE == "FV") then
  set fv_constraints = 1
else
  set fv_constraints = 0
endif

perl $PWD/code/load_balance.pl $TARGET_TASKS $results_dir $PWD $case $casebase $case_dir $fv_constraints $res 

