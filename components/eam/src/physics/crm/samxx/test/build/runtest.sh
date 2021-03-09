#!/bin/bash

ntasks=1
if [[ ! "$1" == "" ]]; then
  ntasks=$1
fi

printf "\nRebuilding\n\n"

make -j8 || exit -1

################################################################################
################################################################################

printf "\n\nRunning 2-D tests\n\n"

printf "\nRunning Fortran code\n\n"
cd fortran2d
rm -f fortran_output_000001.nc
mpirun -n $ntasks ./fortran2d || exit -1
cd ..

printf "\nRunning C++ code\n\n"
cd cpp2d
rm -f cpp_output_000001.nc
mpirun -n $ntasks ./cpp2d || exit -1
cd ..

printf "\nComparing results\n\n"
python nccmp.py fortran2d/fortran_output_000001.nc cpp2d/cpp_output_000001.nc || exit -1

################################################################################
################################################################################

printf "\n\nRunning 3-D tests\n\n"

printf "\nRunning Fortran code\n\n"
cd fortran3d
rm -f fortran_output_000001.nc
mpirun -n $ntasks ./fortran3d || exit -1
cd ..

printf "\nRunning C++ code\n\n"
cd cpp3d
rm -f cpp_output_000001.nc
mpirun -n $ntasks ./cpp3d || exit -1
cd ..

printf "\nComparing results\n\n"
python nccmp.py fortran3d/fortran_output_000001.nc cpp3d/cpp_output_000001.nc || exit -1

################################################################################
################################################################################
