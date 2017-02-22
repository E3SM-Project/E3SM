#!/bin/sh

# Get/Generate the Dashboard Model
if [ $# -eq 0 ]; then
	model=Experimental
else
	model=$1
fi

module purge
module load PrgEnv-intel
module load craype-ivybridge
module load cray-shmem
module load cray-mpich
module load torque
module load git/2.4.6 
module load cmake/3.0.0
module load cray-hdf5-parallel/1.8.14
module load cray-netcdf-hdf5parallel/4.3.3.1
module load cray-parallel-netcdf/1.6.0

export CC=cc
export FC=ftn

export PIO_DASHBOARD_ROOT=`pwd`/dashboard
export PIO_COMPILER_ID=Intel-`$CC --version | head -n 1 | cut -d' ' -f3`

if [ ! -d "$PIO_DASHBOARD_ROOT" ]; then
  mkdir "$PIO_DASHBOARD_ROOT"
fi
cd "$PIO_DASHBOARD_ROOT"

if [ ! -d src ]; then
  git clone https://github.com/PARALLELIO/ParallelIO src
fi
cd src

ctest -S CTestScript.cmake,${model} -VV
