#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps xl netcdf netcdf-fortran cmake python/3.7.0-anaconda3-5.3.0

unset YAKL_ARCH
unset NCRMS

export NCHOME=${OLCF_NETCDF_ROOT}
export NFHOME=${OLCF_NETCDF_FORTRAN_ROOT}
export NCRMS=42
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS="-O3 -qnosave"
export YAKL_CXX_FLAGS=-O3
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"


