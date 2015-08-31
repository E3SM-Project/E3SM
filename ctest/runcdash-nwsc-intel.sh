#!/bin/sh

module reset
module unload netcdf
module swap intel intel/15.0.3
module load git/2.3.0
module load cmake/3.0.2
module load netcdf-mpi/4.3.3.1
module load pnetcdf/1.6.0

export PIO_DASHBOARD_ROOT=`pwd`/dashboard
export PIO_COMPILER_ID=Intel-15.0.3

export CC=mpicc
export FC=mpif90

if [ ! -d "$PIO_DASHBOARD_ROOT" ]; then
  mkdir "$PIO_DASHBOARD_ROOT"
fi
cd "$PIO_DASHBOARD_ROOT"

if [ ! -d src ]; then
  git clone https://github.com/PARALLELIO/ParallelIO src
fi
cd src

ctest -S CTestScript.cmake,Experimental -VV
