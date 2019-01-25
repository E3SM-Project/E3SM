#!/bin/bash

do_cmake=0
do_clean=0
do_make=1
 
HOMME_ROOT=`pwd`/../..
NTRACERS=3
NLEVELS=26
MACH=build_thatchroof

echo $HOMME_ROOT
source ./env_mach_specific
mkdir -p $MACH
cd $MACH

if [ $do_cmake -eq 1 ]; then
NCLIBS="`/opt/netcdf-4.3.3.1_gnu/bin/nf-config --flibs` `/opt/netcdf-4.3.3.1_gnu/bin/nc-config --libs`"
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
  -DNETCDF_DIR="/opt/netcdf-4.3.3.1_gnu"                                         \
  -DPNETCDF_DIR="/opt/parallel-netcdf-1.6.1_gnu"                                 \
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
  -DLINKER_TAG="${NCLIBS}"                             \
  $HOMME_ROOT                                                       || exit -1                                              
fi

if [ $do_clean -eq 1 ]; then
make clean                                                          || exit -1
fi

if [ $do_make -eq 1 ]; then
make -j8 theta-l                                                   || exit -1
make -j8 theta-l-acc                                               || exit -1
fi


