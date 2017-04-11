#!/bin/sh

# Get/Generate the Dashboard Model
if [ $# -eq 0 ]; then
	model=Experimental
else
	model=$1
fi
module rm PrgEnv-intel
module load craype-network-aries
module load PrgEnv-cray
module load craype-haswell
module load cray-shmem
module load cray-mpich
module swap cce cce/8.4.2
module load git/2.6.3
module load cmake/3.3.2
module load cray-hdf5-parallel/1.8.14
module load cray-netcdf-hdf5parallel/4.3.3.1
module load cray-parallel-netcdf/1.6.1

export CC=cc
export FC=ftn

export PIO_DASHBOARD_ROOT=`pwd`/dashboard
export PIO_COMPILER_ID=Cray-`$CC -V 2>&1 | cut -d' ' -f5`

if [ ! -d "$PIO_DASHBOARD_ROOT" ]; then
  mkdir "$PIO_DASHBOARD_ROOT"
fi
cd "$PIO_DASHBOARD_ROOT"

if [ ! -d src ]; then
  git clone https://github.com/PARALLELIO/ParallelIO src
fi
cd src

ctest -S CTestScript.cmake,${model} -VV
