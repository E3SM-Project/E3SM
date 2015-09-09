#!/bin/sh

# Get/Generate the Dashboard Model
if [ $# -eq 0 ]; then
	model=Experimental
else
	model=$1
fi

# Manually set environment variables for CTest run/build
GIT=/soft/versioning/git/2.3.0/bin/git
CTEST=/soft/buildtools/cmake/3.3.0/bin/ctest

export LIBZ=/soft/libraries/alcf/current/xl/ZLIB/
export HDF5=/soft/libraries/hdf5/current/cnk-xl/current
export NETCDF=/soft/libraries/netcdf/current/cnk-xl/current
export PNETCDF=/soft/libraries/pnetcdf/current/cnk-xl/current

export CC=/soft/compilers/wrappers/xl/mpicc
export FC=/soft/compilers/wrappers/xl/mpif90

export PIO_DASHBOARD_ROOT=`pwd`/dashboard
export PIO_COMPILER_ID=Cray-`$CC -qversion | head -n 2 | tail -n 1 | cut -d' ' -f2`

if [ ! -d "$PIO_DASHBOARD_ROOT" ]; then
	mkdir "$PIO_DASHBOARD_ROOT"
fi
cd "$PIO_DASHBOARD_ROOT"

if [ ! -d src ]; then
	$GIT clone https://github.com/PARALLELIO/ParallelIO src
fi
cd src

$CTEST -S CTestScript.cmake,${model} -VV
