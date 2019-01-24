#!/bin/bash

do_cmake=0
do_clean=0
do_make=1
 
HOMME_ROOT=`pwd`/../..
NTRACERS=10
NLEVELS=128
MACH=build_summit

echo $HOMME_ROOT
source ./env_mach_specific
mkdir -p $MACH
cd $MACH

if [ $do_cmake -eq 1 ]; then
NCLIBS="`nf-config --flibs` `nc-config --libs`"
rm -rf CMakeFiles CMakeCache.txt
cmake                                                                            \
  -C $HOMME_ROOT/cmake/machineFiles/darwin.cmake                                 \
  -DCMAKE_Fortran_COMPILER=mpif90                                                \
  -DCMAKE_C_COMPILER=mpicc                                                       \
  -DOPT_FLAGS="-O3"                                                              \
  -DOPENACC_Fortran_FLAGS="-ta=nvidia,cc70,ptxinfo -Minfo=accel"                 \
  -DDEBUG_FLAGS=" "                                                              \
  -DPREQX_NP=4                                                                   \
  -DQSIZE_D=$NTRACERS                                                            \
  -DPREQX_PLEV=$NLEVELS                                                          \
  -DNETCDF_DIR="/ccs/home/imn/netcdf-summit-pgi18.7"                             \
  -DPNETCDF_DIR="$OLCF_PARALLEL_NETCDF_ROOT"                                     \
  -DCMAKE_EXE_LINKER_FLAGS="$NCLIBS"                                             \
  -DBUILD_HOMME_SWEQX=FALSE                                                      \
  -DENABLE_OPENMP=TRUE                                                           \
  -DOPENACC_Linker_FLAGS="-ta=nvidia,cc70"                                       \
  -DBUILD_HOMME_SWEQX=FALSE                                                      \
  -DBUILD_HOMME_PREQX=FALSE                                                      \
  -DBUILD_HOMME_THETA=TRUE                                                       \
  -DBUILD_HOMME_PREQX_ACC=FALSE                                                  \
  -DBUILD_HOMME_PREQX_KOKKOS=FALSE                                               \
  -DBUILD_HOMME_PESE=FALSE                                                       \
  -DBUILD_HOMME_SWIM=FALSE                                                       \
  -DBUILD_HOMME_PRIM=FALSE                                                       \
  $HOMME_ROOT                                                       || exit -1                                              
fi

if [ $do_clean -eq 1 ]; then
make clean                                                          || exit -1
fi

if [ $do_make -eq 1 ]; then
make -j8 theta-l                                                   || exit -1
make -j8 theta-l-acc                                               || exit -1
fi


