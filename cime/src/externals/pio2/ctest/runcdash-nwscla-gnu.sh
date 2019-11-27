#!/bin/sh

# Get/Generate the Dashboard Model
if [ $# -eq 0 ]; then
	model=Experimental
else
	model=$1
fi

module reset
module unload netcdf
module swap intel gnu/8.1.0
module swap mpt openmpi/3.1.0
module load git/2.10.2
module load cmake/3.12.1
module load netcdf/4.6.1
module load pnetcdf/1.10.0

export CC=mpicc
export FC=mpif90

export PIO_DASHBOARD_ROOT=/glade/u/home/jedwards/sandboxes/dashboard
export PIO_COMPILER_ID=GNU-`$CC --version | head -n 1 | tail -n 1 | cut -d' ' -f3`

if [ ! -d "$PIO_DASHBOARD_ROOT" ]; then
  mkdir "$PIO_DASHBOARD_ROOT"
fi
cd "$PIO_DASHBOARD_ROOT"

if [ ! -d src ]; then
  git clone https://github.com/PARALLELIO/ParallelIO src
fi
cd src
git checkout master
git pull origin master

ctest -S CTestScript.cmake,${model} -VV
