#!/bin/sh

# source modules.aurora.oneapi

./configure \
--with-cc=/soft/restricted/CNDA/updates/mpich/52.2/mpich-ofi-all-icc-default-pmix-gpu-drop52/bin/mpicc \
--with-cxx=/soft/restricted/CNDA/updates/mpich/52.2/mpich-ofi-all-icc-default-pmix-gpu-drop52/bin/mpicxx \
--with-fc=/soft/restricted/CNDA/updates/mpich/52.2/mpich-ofi-all-icc-default-pmix-gpu-drop52/bin/mpif90 \
--with-fortran-bindings=1 \
--with-mpiexec="mpiexec " \
--with-batch=0 \
--download-kokkos \
--download-kokkos-commit=origin/develop \
--download-kokkos-kernels-commit=origin/develop \
--download-kokkos-kernels \
--with-kokkos-kernels-tpl=0 \
--with-make-np=8 \
--with-64-bit-indices=0 \
--download-netcdf \
--download-pnetcdf \
--download-hdf5 \
--download-parmetis \
--download-metis \
--download-zlib \
--download-scalapack \
--download-sowing \
--download-triangle \
--download-exodusii \
--download-libceed \
--with-debugging=0 \
--with-sycl=1 \
--with-syclc=/soft/compilers/oneapi/2023.12.15.001/oneapi/compiler/2024.0/bin/icpx \
--with-sycl-arch=pvc \
--SYCLPPFLAGS="-Wno-tautological-constant-compare " \
--COPTFLAGS="-O2 -g " \
--FOPTFLAGS="-O2 -g " \
--CXXOPTFLAGS="-O2 -g " \
--SYCLOPTFLAGS="-O2 -g " \
PETSC_ARCH=aurora-opt-32bit-oneapi-ifx-fc288817

