#!/bin/bash

# Path to ACME and HOMME
export ACME_ROOT=$HOME/workspace/acme
export HOMME_ROOT=$ACME_ROOT/components/homme

# HOMME Settings
NP=4
NLEVELS=26
NTRACERS=4

# machine specific settings (compilers, libraries, etc.)
source /usr/local/tools/dotkit/init.sh

use cmake-3.4.1
use icc-17.0.174
use mvapich2-intel-2.2
use hdf5-intel-parallel-mvapich2-1.10.0
use netcdf-intel-4.3.3.1
use netcdf-fortran-intel-4.4.2

# make build directory
rm -rf $HOMME_ROOT/build
mkdir -p $HOMME_ROOT/build || exit -1
cd $HOMME_ROOT/build

# configure build
cmake \
    \
    -D CMAKE_Fortran_COMPILER=mpif90 \
    -D CMAKE_C_COMPILER=mpicc        \
    -D CMAKE_CXX_COMPILER=mpicxx     \
    \
    -D OPT_FLAGS="-O0"  \
    -D DEBUG_FLAGS="-g" \
    \
    -D NETCDF_DIR=$NETCDF            \
    -D WITH_PNETCDF=FALSE            \
    -D HDF5_DIR=$HDF5                \
    -D HOMME_FIND_BLASLAPACK=TRUE    \
    \
    -D PREQX_NP=$NP                  \
    -D PREQX_PLEV=$NLEVELS           \
    -D QSIZE_D=$NTRACERS             \
    \
    -D BUILD_HOMME_SWEQX=OFF         \
    -D BUILD_HOMME_PREQX=OFF         \
    -D BUILD_HOMME_THETA=ON          \
    -D BUILD_HOMME_PREQX_ACC=OFF     \
    -D BUILD_HOMME_PESE=OFF          \
    -D BUILD_HOMME_SWIM=OFF          \
    -D BUILD_HOMME_PRIM=OFF          \
    \
    -D ENABLE_OPENMP=FALSE           \
    -D ENABLE_OPENACC=FALSE          \
    \
    -D HOMME_PROJID=climalg          \
    \
    $HOMME_ROOT

# build HOMME
make -j $1 theta || exit -1
