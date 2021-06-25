#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps pgi/20.4 netcdf netcdf-fortran cmake python/3.7.0-anaconda3-5.3.0 pgi-cxx14

unset ARCH
unset NCRMS

export PGI_ACC_POOL_ALLOC_MAXSIZE=64B

export NCHOME=${OLCF_NETCDF_ROOT}
export NFHOME=${OLCF_NETCDF_FORTRAN_ROOT}
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS="-O3 -Mextend -ta=nvidia,cc70,ptxinfo,fastmath,managed"
export CXXFLAGS="-O3 -DUSE_ORIG_FFT"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_CUB_HOME="/ccs/home/$USER/cub"

