#!/bin/bash

#export LD_PRELOAD=/usr/lib64/libpil4dfs.so
#printenv
echo "============== MPIEXEC for AURORA ==================="
echo "Launching job after pre-loading the DAOS Interception library..."
#mpiexec $@
mpiexec -genv LD_PRELOAD=/usr/lib64/libpil4dfs.so --no-vni $@
#mpiexec -genv LD_PRELOAD=/soft/perftools/darshan/darshan-3.4.7/lib/libdarshan.so:/opt/aurora/24.347.0/spack/unified/0.9.2/install/linux-sles15-x86_64/oneapi-2025.0.5/hdf5-1.14.5-zrlo32i/lib/libhdf5.so:/opt/aurora/24.347.0/spack/unified/0.9.2/install/linux-sles15-x86_64/oneapi-2025.0.5/parallel-netcdf-1.12.3-cszcp66/lib/libpnetcdf.so:/usr/lib64/libpil4dfs.so --no-vni $@
