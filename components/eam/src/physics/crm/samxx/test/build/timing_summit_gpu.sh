#!/bin/bash

source summit_env
if [[ "$1" != "" ]] ; then
export NCRMS=$1
fi
./cmakescript.sh crmdata_nx32_ny1_nz28_nxrad2_nyrad1_80x.nc crmdata_nx8_ny8_nz28_nxrad2_nyrad2_80x.nc
make -j cpp3d fortran3d || exit -1

printf "Running 3-D tests\n\n"

printf "Running Fortran code\n\n"
cd fortran3d
rm -f fortran_output_000001.nc
jsrun -n 2 -a 21 -c 21 ./fortran3d || exit -1

printf "Running C++ code\n\n"
cd ../cpp3d
rm -f cpp_output_000001.nc
jsrun -n 6 -a 1 -c 1 -g 1 ./cpp3d || exit -1


