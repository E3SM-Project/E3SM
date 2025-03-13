#!/bin/sh

# source modules.compy.intel

./configure \
--with-mpi-dir=/share/apps/intel/2020/compilers_and_libraries_2020.0.166/linux/mpi/intel64 \
--with-netcdf-dir=/share/apps/netcdf/4.6.3/intel/20.0.0 \
--with-pnetcdf-dir=/share/apps/pnetcdf/1.9.0/intel/20.0.0/intelmpi/2020 \
--with-hdf5-dir=/share/apps/hdf5/1.10.5/serial/ \
--with-fortran-bindings=1 \
--download-kokkos \
--download-kokkos-kernels \
--with-kokkos-kernels-tpl=0 \
--with-make-np=8 \
--with-64-bit-indices=1 \
--download-hypre=yes \
--download-mumps=yes \
--download-blas=yes \
--download-lapack=yes \
--download-fblaslapack=yes \
--download-parmetis \
--download-metis \
--download-muparser \
--download-zlib \
--download-scalapack \
--download-sowing \
--download-triangle \
--download-exodusii \
--download-libceed \
--download-cgns-commit=HEAD \
--with-debugging=1 \
PETSC_ARCH=opt-64bit-v3.22.0

