#!/bin/sh

# Get/Generate the Dashboard Model
if [ $# -eq 0 ]; then
	model=Experimental
else
	model=$1
fi

module purge
module load compiler/nag/6.0
module load  tool/parallel-netcdf/1.6.1/nag/openmpi

export CC=mpicc
export FC=mpif90
export PIO_DASHBOARD_SITE="CGD"
export PIO_DASHBOARD_ROOT=/scratch/cluster/jedwards/dashboard
export PIO_COMPILER_ID=Nag-6.0-gcc-`gcc --version | head -n 1 | cut -d' ' -f3`

if [ ! -d "$PIO_DASHBOARD_ROOT" ]; then
  mkdir "$PIO_DASHBOARD_ROOT"
fi
cd "$PIO_DASHBOARD_ROOT"

if [ ! -d src ]; then
  git clone https://github.com/PARALLELIO/ParallelIO src
fi

cd src

ctest -S CTestScript.cmake,${model} -VV
