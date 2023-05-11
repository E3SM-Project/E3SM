#!/bin/bash

unset YAKL_ARCH
unset NCRMS

export NCHOME=${NETCDF_PATH}
export NFHOME=${NETCDF_PATH}
export NCRMS=42
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS="-O3 -ffree-line-length-none -I${NETCDF_PATH}/include"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_ARCH="CUDA"
export YAKL_CUDA_FLAGS="-arch sm_35 -O3 --use_fast_math"


