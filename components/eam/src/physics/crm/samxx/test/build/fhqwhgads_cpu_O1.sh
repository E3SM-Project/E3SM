#!/bin/bash

unset ARCH
unset NCRMS

export NCHOME=/opt/netcdf_gnu
export NFHOME=/opt/netcdf_gnu
export NCRMS=42
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS="-O1 -ffree-line-length-none"
export CXXFLAGS="-O1"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_CUB_HOME="/ccs/home/$USER/cub"


