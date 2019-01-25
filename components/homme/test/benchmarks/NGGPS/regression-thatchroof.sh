#!/bin/bash

build_dir=`pwd`/../../../compile_scripts/thatchroof
test_dir=`pwd`

cd $build_dir
./compile.cmake || exit -1

cd $test_dir
./thatchroof-cpu.job | tee out_cpu.txt
./thatchroof-acc.job | tee out_acc.txt

./nccmp.py ~/homme-runs/cpu/laptop/movies/asp_baroclinic1.nc ~/homme-runs/acc/laptop/movies/asp_baroclinic1.nc | tee diff.txt

