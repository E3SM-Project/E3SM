#!/bin/bash

if [[ "$1" != "" ]] ; then
export NCRMS=$1
fi
./cmakescript.sh crmdata_nx32_ny1_nz28_nxrad2_nyrad1_60x.nc crmdata_nx8_ny8_nz28_nxrad2_nyrad2_60x.nc
make -j
cd fortran2d

printf "Running 2-D tests\n\n"

printf "Running Fortran code\n\n"
cd ../fortran2d
rm -f fortran_output_000001.nc
/gpfs/alpine/world-shared/cli115/mpirun.summit -n 84 -N 84 ./fortran2d

printf "Running C++ code\n\n"
cd ../cpp2d
rm -f cpp_output_000001.nc
jsrun -n 6 -a 1 -c 1 -g 1 ./cpp2d || exit -1
# /gpfs/alpine/world-shared/cli115/mpirun.summit -gpu -n 6 -N 6 ./cpp2d

printf "Running 3-D tests\n\n"

printf "Running Fortran code\n\n"
cd ../fortran3d
rm -f fortran_output_000001.nc
/gpfs/alpine/world-shared/cli115/mpirun.summit -n 84 -N 84 ./fortran3d

printf "Running C++ code\n\n"
cd ../cpp3d
rm -f cpp_output_000001.nc
jsrun -n 6 -a 1 -c 1 -g 1 ./cpp3d || exit -1
# /gpfs/alpine/world-shared/cli115/mpirun.summit -gpu -n 6 -N 6 ./cpp3d


