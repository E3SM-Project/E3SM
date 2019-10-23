#!/bin/sh

### USER OPTION - CHANGE the following to the actual source code location ###
SCREAM_SRC_LOC=~/scream

### Begin Build code
module load gcc
module load cuda
module load cmake
module load spectrum-mpi 
module load netcdf

make -j8 clean
rm -rf CMake*

MPI_F_COMPILE=`mpifort --showme:compile`
MPI_C_COMPILE=`mpicc --showme:compile`
MPI_F_LINK=`mpifort --showme:link`
MPI_C_LINK=`mpicc --showme:link`

cmake \
      -DCMAKE_EXE_LINKER_FLAGS:STRING="${MPI_F_LINK} ${MPI_C_LINK}"                     \
      -DCMAKE_C_FLAGS:STRING="${MPI_C_COMPILE}"                              \
      -DCMAKE_CXX_FLAGS:STRING="${MPI_C_COMPILE} -std=c++11"                            \
      -DCMAKE_Fortran_FLAGS:STRING="${MPI_F_COMPILE}"                                   \
      -D CMAKE_C_COMPILER=mpiCC                                                         \
      -D CMAKE_CXX_COMPILER:STRING=${SCREAM_SRC_LOC}/externals/kokkos/bin/nvcc_wrapper  \
      -D CMAKE_Fortran_COMPILER=mpif90                                                  \
      -D SCREAM_DOUBLE_PRECISION:BOOL=ON                                                \
      -D KOKKOS_ENABLE_DEBUG:BOOL=FALSE                                                 \
      -D KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION=FALSE                                   \
      -D KOKKOS_ENABLE_CUDA_LAMBDA=TRUE                                                 \
      -D KOKKOS_ENABLE_DEPRECATED_CODE=FALSE                                            \
      -D KOKKOS_ENABLE_EXPLICIT_INSTANTIATION=FALSE                                     \
      -D KOKKOS_ENABLE_OPENMP:BOOL=Off                                                  \
      -D KOKKOS_ENABLE_CUDA:BOOL=ON                                                     \
      -D KOKKOS_ARCH:STRING="Volta70"                                                   \
      ${SCREAM_SRC_LOC}/components/scream

VERBOSE=1 make -j8
