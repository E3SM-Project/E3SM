#!/bin/bash

source $MODULESHOME/init/bash
module load PrgEnv-amd craype-accel-amd-gfx90a rocm cray-hdf5 cray-netcdf

unset YAKL_ARCH
unset NCRMS
unset CXXFLAGS
unset FFLAGS
unset FCLAGS

export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_SMP_SINGLE_COPY_MODE=CMA

export NCHOME=${NETCDF_DIR}
export NFHOME=/ccs/home/imn/netcdf-fortran-crusher-gnu
export NCRMS=42
export CC=hipcc
export CXX=hipcc
export FC=gfortran
export YAKL_F90_FLAGS="-O3 -ffree-line-length-none"
export FFLAGS="-O3 -ffree-line-length-none -I${NFHOME}/include"
export YAKL_C_FLAGS="-O3"
export YAKL_HIP_FLAGS="-O3 -Wno-tautological-pointer-compare -Wno-unused-result -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx90a -x hip"
export YAKL_ARCH="HIP"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"


