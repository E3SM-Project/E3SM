#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc cuda/11.4.0 netcdf-c netcdf-fortran cmake python/3.7.0-anaconda3-5.3.0

unset YAKL_ARCH
unset NCRMS
unset CXXFLAGS
unset FFLAGS
unset FCLAGS

export NCHOME=${OLCF_NETCDF_C_ROOT}
export NFHOME=${OLCF_NETCDF_FORTRAN_ROOT}
export NCRMS=42
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export YAKL_F90_FLAGS="-O3 -ffree-line-length-none"
export FFLAGS="-O3 -ffree-line-length-none -I${OLCF_NETCDF_FORTRAN_ROOT}/include"
export YAKL_C_FLAGS="-O3"
export YAKL_CUDA_FLAGS="-arch sm_70 -O3 --use_fast_math"
export YAKL_ARCH="CUDA"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"


