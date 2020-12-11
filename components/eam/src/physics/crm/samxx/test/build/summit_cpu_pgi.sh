#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps pgi netcdf netcdf-fortran cmake python/3.7.0-anaconda3-5.3.0

unset ARCH
unset NCRMS

export NCHOME=${OLCF_NETCDF_ROOT}
export NFHOME=${OLCF_NETCDF_FORTRAN_ROOT}
export NCRMS=84
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS="-O3 -Mextend"
export CXXFLAGS="-O3 -DUSE_ORIG_FFT"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_CUB_HOME="/ccs/home/$USER/cub"

