#!/bin/sh

#source modules.pm-gpu.gnugpu

./configure \
--with-cc=/opt/cray/pe/craype/2.7.31.11/bin/cc \
--with-cxx=/opt/cray/pe/craype/2.7.31.11/bin/CC \
--with-fc=/opt/cray/pe/craype/2.7.31.11/bin/ftn \
--with-fortran-bindings=1 \
--with-mpiexec="srun -g 8 --smpiargs=-gpu " \
--with-batch=0 \
--with-make-np=8 \
--with-netcdf-dir=/opt/cray/pe/netcdf-hdf5parallel/4.9.0.13/gnu/12.3/ \
--with-pnetcdf-dir=/opt/cray/pe/parallel-netcdf/1.12.3.9/gnu/12.3 \
--with-hdf5-dir=/opt/cray/pe/hdf5-parallel/1.14.3.1/gnu/12.3/ \
--with-hip=1 \
--with-hipc=/opt/rocm-5.4.0/bin/hipcc \
--download-parmetis \
--download-metis \
--download-muparser \
--download-zlib \
--download-scalapack \
--download-sowing \
--download-triangle \
--download-exodusii \
--download-libceed \
--with-debugging=1 \
--with-64-bit-indices=0 \
PETSC_ARCH=frontier-gpu-debug-32bit-gcc-12-3-0-0d6defa7a01

