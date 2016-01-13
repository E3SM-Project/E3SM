#!/usr/bin/env bash

# After this executes, do:
#   make -j 8

# Set path to top cism directory
# Note, this is an easy way to build out of source.
# In directory you want to build in, run:
#   $ source $CISM/builds/linux-gnu-cism/linux-gnu-cism-cmake $CISM
# where $CISM is the path to the top level cism directory.
if [ $# -eq 0 ]
then
    cism_top="../.." 
else
    cism_top=${1}
fi

echo CISM: ${cism_top}


export CRAY_CPU_TARGET=istanbul

source $MODULESHOME/init/bash

module unload cmake
module unload cray-hdf5 
module unload cray-hdf5-parallel
module unload netcdf
module unload python
module unload cray-shmem
module unload cray-mpich cray-mpich2
module unload netcdf-hdf5parallel cray-netcdf-hdf5parallel cray-parallel-netcdf cray-hdf5-parallel
module unload boost 
module unload pgi
module unload PrgEnv-cray PrgEnv-gnu PrgEnv-intel PrgEnv-pathscale PrgEnv-pgi

module load modules
module load cmake/2.8.11.2
module load PrgEnv-pgi/5.2.40
module unload pgi
module load cray-shmem
module load cray-mpich
module load pgi/14.2.0
module load cray-hdf5-parallel/1.8.13
module load cray-netcdf-hdf5parallel/4.3.2
module load python
module load boost/1.57.0

echo module list

# remove old build data:
rm -f ./CMakeCache.txt
rm -rf ./CMakeFiles

# BUILD OPTIONS:
# The call to cmake below includes several input ON/OFF switch parameters, to
# provide a simple way to select different build options.  These are:
# CISM_BUILD_CISM_DRIVER -- ON by default, set to OFF to only build the CISM libraries.
# CISM_ENABLE_BISICLES -- OFF by default, set to ON to build a BISICLES-capable cism_driver.
# CISM_ENABLE_FELIX -- OFF by default, set to ON to build a FELIX-capable cism_driver.
# CISM_USE_TRILINOS -- OFF by default, set to on for builds with Trilinos.
# CISM_MPI_MODE -- ON by default, only set to OFF for serial builds.
# CISM_SERIAL_MODE -- OFF by default, set to ON for serial builds.
# CISM_USE_GPTL_INSTRUMENTATION -- ON by default, set to OFF to not use GPTL instrumentation.
# CISM_COUPLED -- OFF by default, set to ON to build with CESM.

cmake \
  -D CISM_BUILD_CISM_DRIVER:BOOL=ON \
  -D CISM_ENABLE_BISICLES=OFF \
  -D CISM_ENABLE_FELIX=OFF \
\
  -D CISM_USE_TRILINOS:BOOL=${CISM_USE_TRILINOS:=ON} \
  -D CISM_MPI_MODE:BOOL=ON \
  -D CISM_SERIAL_MODE:BOOL=OFF \
\
  -D CISM_USE_GPTL_INSTRUMENTATION:BOOL=ON \
  -D CISM_COUPLED:BOOL=OFF \
\
  -D CISM_TRILINOS_DIR=/lustre/atlas/world-shared/cli900/cesm/software/Trilinos/Trilinos-11.12.1_gptl/titan-pgi-ci-nophal/install \
  -D CISM_TRILINOS_GPTL_DIR=/lustre/atlas/world-shared/cli900/cesm/software/Trilinos/Trilinos-11.12.1_gptl/titan-pgi-ci-nophal/install \
  -D CISM_TRILINOS_ALBANY_DIR=/lustre/atlas/world-shared/cli900/cesm/software/Trilinos/Trilinos-11.12.1_gptl/titan-pgi-ci-nophal/install \
\
  -D CISM_GPTL_DIR=/lustre/atlas/world-shared/cli900/cesm/software/libgptl/libgptl-titan-pgi \
  -D CISM_NETCDF_DIR=/opt/cray/netcdf-hdf5parallel/4.3.2/PGI/141 \
\
  -D CMAKE_INSTALL_PREFIX:PATH=$cism_top/builds/titan-gnu/install \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
  -D CMAKE_VERBOSE_CONFIGURE:BOOL=ON \
\
  -D CMAKE_CXX_COMPILER=CC \
  -D CMAKE_C_COMPILER=cc \
  -D CMAKE_Fortran_COMPILER=ftn \
\
  -D CMAKE_CXX_FLAGS:STRING="-fast -Kieee --diag_suppress 554,111,611" \
  -D CISM_Fortran_FLAGS:STRING="-fast -Kieee" \
  -D CISM_FMAIN=/opt/pgi/14.2.0/linux86-64/14.2/lib/f90main.o \
  -D BISICLES_LIB_SUBDIR=libpgi \
  -D BISICLES_INTERFACE_DIR=$cism_top/../BISICLES/CISM-interface/interface \
  -D CISM_MPI_LIBS:STRING="mpichf90" \
  -D CISM_USE_CXX_IMPLICIT_LIBS:BOOL=OFF \
  -D CISM_STATIC_LINKING:BOOL=ON \
  ${cism_top}

