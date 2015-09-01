#!/bin/sh

# Get/Generate the Dashboard Model
if [ $# -eq 0 ]; then
	model=Experimental
else
	model=$1
fi

# How to "resoft" a specific environment for the dashboard run?
exit 1
#------------------------------------------------------------------------------
export CC=mpicc
export FC=mpif90

export PIO_DASHBOARD_ROOT=`pwd`/dashboard
export PIO_COMPILER_ID=Cray-`$CC -qversion | head -n 2 | tail -n 1 | cut -d' ' -f2`

if [ ! -d "$PIO_DASHBOARD_ROOT" ]; then
  mkdir "$PIO_DASHBOARD_ROOT"
fi
cd "$PIO_DASHBOARD_ROOT"

if [ ! -d src ]; then
  git clone https://github.com/PARALLELIO/ParallelIO src
fi
cd src

ctest -S CTestScript.cmake,${model} -VV
