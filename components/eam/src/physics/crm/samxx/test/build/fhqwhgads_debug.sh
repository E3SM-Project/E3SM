#!/bin/bash

unset ARCH
unset NCRMS

export NCHOME=/opt/netcdf_gnu
export NFHOME=/opt/netcdf_gnu
export NCRMS=42
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS="-O0 -fcheck=all -ffree-line-length-none -g"
export CXXFLAGS="-O0 -DYAKL_DEBUG -g"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_CUB_HOME="/ccs/home/$USER/cub"


