#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load rocm/3.5.0 cmake openmpi

./cmakeclean.sh

export PATH=/ccs/home/imn/libs_gnu8.1.0/netcdf-4.3.3.1_gnu/bin:${PATH}

unset GATOR_DISABLE
unset ARCH
unset NCRMS

export OMPI_CXX=hipcc
export OMPI_CC=hipcc
export NCHOME=/ccs/home/imn/libs_gnu8.1.0/netcdf-4.3.3.1_gnu
export NFHOME=/ccs/home/imn/libs_gnu8.1.0/netcdf-4.3.3.1_gnu
export NCRMS=42
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS="-O3 -ffree-line-length-none"
export CXXFLAGS="-O3"
export HIPFLAGS="-O3"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"

