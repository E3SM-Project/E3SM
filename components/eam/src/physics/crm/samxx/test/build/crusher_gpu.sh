#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load PrgEnv-cray/8.3.3 craype-accel-amd-gfx90a rocm/5.1.0
module load craype-network-ofi cce/14.0.0  cray-mpich/8.1.16
module load cray-hdf5 cray-netcdf cmake

unset ARCH
unset YAKL_ARCH
unset NCRMS
unset MACH
unset CC
unset CXX
unset FC
unset YAKL_HIP_FLAGS
unset YAKL_DEBUG
unset FFLAGS
unset CXXFLAGS

export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_SMP_SINGLE_COPY_MODE=CMA

export MACH="crusher"
export NCHOME=${NETCDF_DIR}
export NFHOME=${NETCDF_DIR}
export MPIHOME=${CRAY_MPICH_DIR}
export NCRMS=224
export CC=hipcc
export CXX=hipcc
export FC=ftn
export FFLAGS=" -O3 -hnoomp -hnoacc "
export CXXFLAGS=" -O3 -I${ROCM_PATH}/include "
export ARCH="HIP"
export YAKL_ARCH="HIP"
export YAKL_HIP_FLAGS="-O3 -munsafe-fp-atomics -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx90a -x hip " 
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_CUB_HOME="`pwd`/../../../../../../../../externals/hipCUB"
