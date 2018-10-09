#!/bin/bash

# Path to HOMME
export HOMME_ROOT=$PWD

# Set build type ("new", "test", or "update" with default being "update")
if [[ $# -gt 4 ]]; then
  BUILD=$5
else
  BUILD="update"
fi

# machine specific settings (compilers, libraries, etc.)
if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
  source /usr/local/tools/dotkit/init.sh
  use cmake-3.4.1
  use icc-17.0.174
  use mvapich2-intel-2.2
  use hdf5-intel-parallel-mvapich2-1.10.0
  use netcdf-intel-4.3.3.1
  use netcdf-fortran-intel-4.4.2
  export SUNDIALS=$HOME/local/sundials-temp_cab_intel_opt
elif [[ $HOSTNAME == "quartz"* ]]; then
  SYSTEM=quartz
  module purge
  module load cmake
  if [[ $BUILD == "test" ]]; then
    module load gcc/6.1.0
    module load openmpi/2.0.0
    export NETCDF=$HOME/local/netcdf-c4.3.3.1_f4.4.2_quartz_gnu_opt
    export HDF5=$HOME/local/hdf5-1.10.1_quartz_gnu_opt
    export SUNDIALS=$HOME/local/sundials-2.7.0_quartz_gnu_opt
  else
    module load intel/16.0.3
    export NETCDF=$HOME/local/netcdf-c-serial-4.5.0_f-4.4.4_quartz_intel_opt
    export HDF5=$HOME/local/hdf5-serial-1.10.3_quartz_intel_opt
    export SUNDIALS=$HOME/local/sundials-3.1.2_quartz_intel_opt
    MPI=$HOME/local/mpich-3.2_quartz_intel_opt/bin
    export MPIF90=$MPI/mpif90
    export MPICC=$MPI/mpicc
    export MPICXX=$MPI/mpicxx
  fi
fi


# HOMME Settings
if [[ $# -gt 1 ]]; then
  NP=$2
  NLEVELS=$3
  NTRACERS=$4
else
  NP=4
  NLEVELS=26
  NTRACERS=4
fi

# Delete and create build directory if necessary
if [[ $BUILD != "update" ]]; then
  rm -rf $HOMME_ROOT/build_${SYSTEM}
fi
mkdir -p $HOMME_ROOT/build_${SYSTEM} || exit -1
cd $HOMME_ROOT/build_${SYSTEM}

# configure build (if necessary)
if [[ $BUILD != "update" ]]; then
  cmake \
    -C $HOMME_ROOT/cmake/machineFiles/lc.cmake \
    \
    -D PREQX_NP=$NP                  \
    -D PREQX_PLEV=$NLEVELS           \
    -D QSIZE_D=$NTRACERS             \
    \
    -D BUILD_HOMME_SWEQX=OFF \
    -D BUILD_HOMME_PREQX=OFF \
    -D BUILD_HOMME_THETA=ON \
    -D BUILD_HOMME_PREQX_ACC=OFF \
    -D BUILD_HOMME_PESE=OFF \
    -D BUILD_HOMME_SWIM=OFF \
    -D BUILD_HOMME_PRIM=OFF\
    \
    -D HOMME_USE_ARKODE=ON \
    \
    $HOMME_ROOT
fi

# build HOMME
make -j$1 theta-l || exit -1
