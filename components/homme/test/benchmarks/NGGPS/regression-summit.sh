#!/bin/bash

build_dir=`pwd`/../../../compile_scripts/summit
test_dir=`pwd`

cd $build_dir
./compile.cmake || exit -1

cd $test_dir
./summit-cpu.job | tee out_cpu.txt
./summit-acc.job | tee out_acc.txt

./nccmp.py /gpfs/alpine/proj-shared/cli115/imn/homme-runs/cpu/tiny/movies/asp_baroclinic1.nc /gpfs/alpine/proj-shared/cli115/imn/homme-runs/acc/tiny/movies/asp_baroclinic1.nc | tee diff.txt

