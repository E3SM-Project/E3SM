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


module unload cmake
module unload hdf5
module unload hdf5-parallel cray-hdf5-parallel
module unload netcdf cray-netcdf-hdf5parallel
module unload python
module unload cray-shmem
module unload cray-mpich2
module unload cray-mpich
module unload boost
module unload gcc
module unload PrgEnv-cray PrgEnv-gnu PrgEnv-intel PrgEnv-pathscale PrgEnv-pgi
module unload craype

module load modules/3.2.10.2
module load cmake/3.1.3
module load PrgEnv-gnu/5.2.40
module load gcc/4.9.2
module load cray-shmem/7.1.1
module load cray-mpich/7.1.1
module load cray-netcdf-hdf5parallel/4.3.2
module load cray-parallel-netcdf/1.5.0
module load python/2.7.9
module load boost/1.57

module list


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
  -D ALBANY_FELIX_DYCORE:BOOL=OFF \
\
  -D CISM_TRILINOS_DIR=/project/projectdirs/piscees/trilinos/trilinos-albany-build/install \
  -D CISM_TRILINOS_GPTL_DIR=/project/projectdirs/piscees/cism_gptl/Trilinos-11.12.1/hopper-gnu-ci-nophal/install \
  -D CISM_TRILINOS_ALBANY_DIR=/project/projectdirs/piscees/trilinos/trilinos-albany-build/install \
\
  -D CISM_NETCDF_DIR=$NETCDF_DIR \
  -D CISM_HDF5_LIB_DIR=$HDF5_DIR \
  -D CISM_MPI_BASE_DIR=$MPICH_DIR \
\
  -D CISM_GPTL_DIR=/project/projectdirs/piscees/cism_gptl/libgptl/libgptl-hopper-gnu \
\
  -D CMAKE_INSTALL_PREFIX:PATH=$cism_top/builds/hopper-gnu/install \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
  -D CMAKE_VERBOSE_CONFIGURE:BOOL=ON \
\
  -D CMAKE_CXX_COMPILER=CC \
  -D CMAKE_C_COMPILER=cc \
  -D CMAKE_Fortran_COMPILER=ftn \
\
  -D CMAKE_CXX_FLAGS:STRING="-O2 -fopenmp" \
  -D CISM_Fortran_FLAGS:STRING="-O2 -fopenmp -ffree-line-length-none " \
  -D BISICLES_LIB_SUBDIR=libgnu \
  -D BISICLES_INTERFACE_DIR=$cism_top/../BISICLES/CISM-interface/interface \
  -D CISM_MPI_LIBS:STRING="mpichf90" \
  -D CISM_STATIC_LINKING:BOOL=ON \
  ${cism_top}

