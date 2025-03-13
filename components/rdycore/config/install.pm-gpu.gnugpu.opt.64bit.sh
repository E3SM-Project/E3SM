#!/bin/sh

# source modules.pm-gpu.gnugpu

./configure \
--with-cc=cc \
--with-cxx=CC \
--with-fc=ftn \
--CFLAGS=" -g " \
--CXXFLAGS=" -g " \
--CUDAFLAGS=" -g -Xcompiler -rdynamic " \
--with-fortran-bindings=1 \
--COPTFLAGS="   -O" \
--CXXOPTFLAGS=" -O" \
--FOPTFLAGS="   -O" \
--with-mpiexec="srun -G4" \
--with-batch=0 \
--download-kokkos \
--download-kokkos-kernels \
--with-kokkos-kernels-tpl=0 \
--with-make-np=8 \
--with-64-bit-indices=1 \
--with-netcdf-dir=/opt/cray/pe/netcdf-hdf5parallel/4.9.0.3/gnu/9.1 \
--with-pnetcdf-dir=/opt/cray/pe/parallel-netcdf/1.12.3.3/gnu/9.1 \
--with-hdf5-dir=/opt/cray/pe/hdf5-parallel/1.12.2.3/gnu/9.1 \
--with-cuda-dir=/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7 \
--with-cuda-arch=80 \
--download-parmetis \
--download-metis \
--download-zlib \
--download-scalapack \
--download-sowing \
--download-triangle \
--download-exodusii \
--download-libceed \
--download-cgns-commit=HEAD \
--with-debugging=0 \
PETSC_ARCH=pm-gpu-opt-64bit-gcc-11-2-0-fc2888174f5

