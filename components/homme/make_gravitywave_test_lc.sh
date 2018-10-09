#!/bin/bash

WORKDIR=$S/homme/gravitywave_test_noheating
NLEV=20

if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
elif [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
else
  echo "Script not implemented for this system"
  exit -1
fi

./cmake_homme_lc.sh 8 4 $NLEV 0 $1

if [[ $? == 0 ]]; then
  mkdir -p $WORKDIR
  rm -f $WORKDIR/theta-l-nlev$NLEV $WORKDIR/submit_gravitywave.py
  rm -f $WORKDIR/run_convergence_test.sh $WORKDIR/dictionaries.py
  rm -f  $WORKDIR/plot_convergence_test.py
  rm -f $WORKDIR/plot_NLsolver_performance.py $WORKDIR/plot_RMSE_efficiency.py
  cp build_$SYSTEM/src/theta-l/theta-l $WORKDIR/theta-l-nlev$NLEV
  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta-l/submit_gravitywave_lc.py $WORKDIR/submit_gravitywave.py
  ln -s $PWD/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta-l/run_convergence_test_lc.sh $WORKDIR/run_convergence_test.sh
  ln -s $PWD/plot_scripts/dictionaries.py $WORKDIR/dictionaries.py
  ln -s $PWD/plot_scripts/plot_convergence_test.py $WORKDIR/plot_convergence_test.py
  ln -s $PWD/plot_scripts/plot_NLsolver_performance.py $WORKDIR/plot_NLsolver_performance.py
  ln -s $PWD/plot_scripts/plot_RMSE_efficiency.py $WORKDIR/plot_RMSE_efficiency.py
fi
