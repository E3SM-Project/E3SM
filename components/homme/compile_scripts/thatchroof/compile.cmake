#!/bin/bash

do_cmake=1
do_clean=1
do_make=1
 
HOMME_ROOT=/home/$USER/ACME/components/homme
NTRACERS=40
NLEVELS=72

source ./env_mach_specific

mkdir -p thatchroof-openacc     || exit -1
cd thatchroof-openacc

NCLIBS="`/opt/netcdf-4.3.3.1_pgi16.4/bin/nf-config --flibs` `/opt/netcdf-4.3.3.1_pgi16.4/bin/nc-config --libs`"

if [ $do_cmake -eq 1 ]; then
rm -rf CMakeFiles CMakeCache.txt
cmake                                                                            \
  -C $HOMME_ROOT/cmake/machineFiles/darwin.cmake                                 \
  -DCMAKE_Fortran_COMPILER=mpif90                                                \
  -DCMAKE_C_COMPILER=mpicc                                                       \
  -DOPT_FLAGS="-O3"                                                              \
  -DOPENACC_Fortran_FLAGS="-ta=nvidia,cuda8.0,cc35,ptxinfo,pinned -Minfo=accel"  \
  -DDEBUG_FLAGS=" "                                                              \
  -DNETCDF_DIR=$NETCDF_DIR                                                       \
  -DPREQX_NP=4                                                                   \
  -DQSIZE_D=$NTRACERS                                                            \
  -DPREQX_PLEV=$NLEVELS                                                          \
  -DBUILD_HOMME_SWEQX=FALSE                                                      \
  -DBUILD_HOMME_PREQX=TRUE                                                       \
  -DBUILD_HOMME_PREQX_ACC=TRUE                                                   \
  -DENABLE_OPENMP=TRUE                                                           \
  -DCMAKE_EXE_LINKER_FLAGS="${NCLIBS}"                                           \
  -DOPENACC_Linker_FLAGS="-ta=nvidia,cuda8.0,cc35,pinned"                        \
  -DAPPEND_LIBRARIES="netcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz -ldl"           \
  $HOMME_ROOT                                                       || exit -1                                              
fi
#APPEND_LIBRARIES is a hack because CMake is really really stupid

if [ $do_clean -eq 1 ]; then
make clean                                                          || exit -1
fi

if [ $do_make -eq 1 ]; then
make -j24 preqx_acc preqx                                           || exit -1
mkdir -p $HOMME_ROOT/build/preqx
cp ./src/preqx_acc/preqx_acc $HOMME_ROOT/build/preqx/preqx.openacc  || exit -1
cp ./src/preqx/preqx         $HOMME_ROOT/build/preqx/preqx          || exit -1
fi
