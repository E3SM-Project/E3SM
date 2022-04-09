#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc/9.1.0 cuda/11.0.3 netcdf-c/4.8.0 netcdf-fortran/4.4.5 cmake  python/3.8-anaconda3 forge/20.1 

unset ARCH
unset NCRMS
unset MACH
unset CC
unset CXX
unset FC

export MACH="summit"
export NCHOME=${OLCF_NETCDF_C_ROOT}
export NFHOME=${OLCF_NETCDF_FORTRAN_ROOT}
export NCRMS=168
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS=" -O3 -ffree-line-length-none "
export CXXFLAGS=" -O3 -DTHRUST_IGNORE_CUB_VERSION_CHECK "
export ARCH="CUDA"
export YAKL_ARCH="CUDA"
export YAKL_CUDA_FLAGS="-arch sm_70 --use_fast_math -D__USE_CUDA__ -DTHRUST_IGNORE_CUB_VERSION_CHECK --expt-extended-lambda --expt-relaxed-constexpr"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_CUB_HOME="`pwd`/../../../../../../../../externals/cub"


