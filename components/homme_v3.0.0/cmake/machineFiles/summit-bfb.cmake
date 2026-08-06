#interactive job
#bsub -W 2:00 -nnodes 1 -P cli115 -Is /bin/bash

#1 gpu on 1 node
#jsrun -n 1 -r 1 -l gpu-gpu -b packed:1 -d plane:1 -a1 -c7 -g1 --smpiargs "-gpu" EXEC  < ${nlname}
#6 gpus on 1 node
#jsrun -n 6 -r 6 -l gpu-gpu -b packed:1 -d plane:1 -a1 -c7 -g1 --smpiargs "-gpu" EXEC  < ${nlname}
#6 ranks, cpu only on 1 node, not 42 rank, NEEDS visible GPU due to some kokkos init, so, -g1, not -g0
#jsrun -n 6 -r 6 -l cpu-cpu -b packed:1 -d plane:1 -a1 -c7 -g1  EXEC  < ${nlname}


set(CMAKE_C_FLAGS "-w" CACHE STRING "")
set(ADD_CXX_FLAGS "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -Wno-unknown-pragmas --fmad=false -O0" CACHE STRING "")
set(ADD_Fortran_FLAGS " -ffp-contract=off -O0" CACHE STRING "")
set(OPT_FLAGS "-O0" CACHE STRING "")
set(DEBUG_FLAGS "-ffp-contract=off -g" CACHE STRING "")


#SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")
SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM "8" CACHE STRING  "")

SET (NETCDF_DIR $ENV{OLCF_NETCDF_FORTRAN_ROOT} CACHE FILEPATH "")
SET (NetCDF_Fortran_PATH $ENV{OLCF_NETCDF_FORTRAN_ROOT}  CACHE STRING "")
SET (NetCDF_C_PATH $ENV{OLCF_NETCDF_C_ROOT}  CACHE STRING "")
#SET(NetCDF_C_LIBRARY "/sw/summit/spack-envs/base/opt/linux-rhel8-ppc64le/gcc-7.5.0/netcdf-c-4.8.0-pwi4jbrnwv4lrrjxdu5czbos5uvvjgvr/lib" CACHE STRING "")
#SET(NetCDF_C_INCLUDE_DIR "/sw/summit/spack-envs/base/opt/linux-rhel8-ppc64le/gcc-7.5.0/netcdf-c-4.8.0-pwi4jbrnwv4lrrjxdu5czbos5uvvjgvr/include" CACHE STRING "")

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE TRUE CACHE BOOL "")
#SET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")

SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

SET(USE_TRILINOS OFF CACHE BOOL "")

SET(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA ON CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
SET(Kokkos_ARCH_VOLTA70 ON CACHE BOOL "")
SET(Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/../../externals/kokkos/bin/nvcc_wrapper" CACHE STRING "")

set (ENABLE_OPENMP OFF CACHE BOOL "")
set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

set (HOMME_TESTING_PROFILE "dev" CACHE STRING "")

set (USE_NUM_PROCS 4 CACHE STRING "")

#set (OPT_FLAGS "-mcpu=power9 -mtune=power9" CACHE STRING "")

set (USE_MPI_RUN_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/cmake/machineFiles/summit-run.sh" CACHE STRING "")
