#!/bin/bash

do_cmake=1
do_clean=1
do_make=1


# Change next two lines to your paths
export HOMME_ROOT=$PROJWORK/cli115/4ue/E3SM/components/homme

source $HOMME_ROOT/compile_scripts/titan/env_mach_specific.cpu.gnu
export Trilinos_DIR=/lustre/atlas1/cli106/proj-shared/4ue/Trilinos/reg_build/install
#export Trilinos_DIR=/opt/cray/trilinos/12.12.1.0/GNU/7.1/x86_64

# leaving this in causes dlopen link errors, leaving it out causes "undefined reference to `main' "
export CRAYPE_LINK_TYPE='dynamic' 

export TARGET='titan-swim-cpu-gnu'

mkdir -p $TARGET     || exit -1
cd $TARGET

if [ $do_cmake -eq 1 ]; then
rm -rf CMakeFiles CMakeCache.txt
cmake                                                                          \
  -C $HOMME_ROOT/cmake/machineFiles/titan.cmake.gnu                            \
  -DCMAKE_BUILD_TYPE=RELEASE                                                   \
  -DCMAKE_Fortran_COMPILER=ftn                                                 \
  -DCMAKE_C_COMPILER=cc                                                        \
  -DCMAKE_CXX_COMPILER=CC                                                      \
  -DOPT_FLAGS="-O2"           			                               \
  -DDEBUG_FLAGS=" "                                                            \
  -DNETCDF_DIR=$NETCDF_DIR                                                     \
  -DWITH_PNETCDF=FALSE                                                         \
  -DHDF5_DIR=$HDF5_DIR                                                         \
  -DSWIM_NP=4                                                                  \
  -DBUILD_HOMME_SWIM=TRUE                                                      \
  -DBUILD_HOMME_SWEQX=FALSE                                                    \
  -DBUILD_HOMME_PRIM=FALSE                                                     \
  -DBUILD_HOMME_PREQX=FALSE                                                    \
  -DENABLE_OPENMP=TRUE                                                         \
  -DHOMME_PROJID=cli106ms                                                      \
  -DENABLE_OPENACC=FALSE                                                       \
  -DENABLE_CUDA_FORTRAN=FALSE                                                  \
  -DHOMME_FIND_BLASLAPACK=TRUE                                                 \
  -DNetcdf_NC_CONFIG_BIN="/opt/cray/netcdf/4.3.3.1/bin"                        \
  $HOMME_ROOT
fi

if [ $do_clean -eq 1 ]; then
make clean                                                  || exit -1
rm -rf $TARGET
fi

if [ $do_make -eq 1 ]; then
#VERBOSE=1 
make -j16 swim                                            || exit -1
mkdir -p $HOMME_ROOT/build/swim
cp ./src/sweqx/swim $HOMME_ROOT/build/swim/swim.cpu.gnu  || exit -1
fi
