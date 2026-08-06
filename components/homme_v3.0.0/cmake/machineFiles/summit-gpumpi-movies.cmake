#interactive job
#bsub -W 2:00 -nnodes 1 -P cli115 -Is /bin/bash


#SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")
SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM "16" CACHE STRING  "")

SET (NETCDF_DIR $ENV{OLCF_NETCDF_FORTRAN_ROOT} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{OLCF_HDF5_ROOT} CACHE FILEPATH "")

#for scorpio
SET (NetCDF_C_PATH $ENV{OLCF_NETCDF_ROOT} CACHE FILEPATH "")

#SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

SET(BUILD_HOMME_SWEQX FALSE CACHE BOOL "")
SET(BUILD_HOMME_PREQX_ACC FALSE CACHE BOOL "")
SET(BUILD_HOMME_PREQX FALSE CACHE BOOL "")
SET(BUILD_HOMME_PREQX_KOKKOS FALSE CACHE BOOL "")
SET(BUILD_HOMME_THETA FALSE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")
SET(BUILD_HOMME_TOOL FALSE CACHE BOOL "")

#SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

SET(USE_TRILINOS OFF CACHE BOOL "")

SET(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA ON CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
SET(Kokkos_ARCH_VOLTA70 ON CACHE BOOL "")
SET(Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA_ARCH_LINKING OFF CACHE BOOL "")

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "/ccs/home/onguba/kokkos/bin/nvcc_wrapper" CACHE STRING "")

set (ENABLE_OPENMP OFF CACHE BOOL "")
set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

set (HOMME_TESTING_PROFILE "dev" CACHE STRING "")

#set (OPT_FLAGS "-mcpu=power9 -mtune=power9" CACHE STRING "")

#temp change to have bngry exchange file compile with cuda
#will get rid of this later
#ADD_Fortran_FLAGS, ADD_C_FLAGS, ADD_CXX_FLAGS
set (ADD_CXX_FLAGS "--save-temps" CACHE STRING "")

SET (USE_MPI_OPTIONS "--bind-to core" CACHE FILEPATH "")


###################################################################

set(HOMME_MACHINE "summit-gpu" CACHE STRING "")

set(SUMMIT_NODES "2" CACHE STRING "")
set(SUMMIT_RES_PER_NODE "6" CACHE STRING "")
math(EXPR SUMMIT_NRES "${SUMMIT_NODES} * ${SUMMIT_RES_PER_NODE}")

set(varnames  CACHE STRING "")
set(varvals  CACHE STRING "")

set(varnames ${varnames}  "SUMMIT_NODES" CACHE INTERNAL "")
set(varvals ${varvals} "${SUMMIT_NODES}" CACHE INTERNAL "")

set(varnames ${varnames}  "SUMMIT_MODULES_GPU" CACHE INTERNAL "")
set(varvals ${varvals} "load-modules-gpu" CACHE INTERNAL "")

set(varnames ${varnames}  "SUMMIT_MODULES_P9" CACHE INTERNAL "")
set(varvals ${varvals} "load-modules-p9" CACHE INTERNAL "")

set(varnames ${varnames}  "SUMMIT_JSRUN_GPU"  CACHE INTERNAL "")
set(varvals ${varvals}  "jsrun -n ${SUMMIT_NRES} -r ${SUMMIT_RES_PER_NODE} -l gpu-gpu -b packed:1 -d plane:1 -a1 -c7 -g1 --smpiargs -gpu " CACHE INTERNAL "")

set(varnames ${varnames}  "SUMMIT_JSRUN_P9" CACHE INTERNAL "")
set(varvals ${varvals}  "jsrun -n ${SUMMIT_NRES} -r ${SUMMIT_RES_PER_NODE} -l cpu-cpu -b packed:1 -d plane:1 -a7 -c7 -g0 "  CACHE INTERNAL "")

set(varnames ${varnames}  "SUMMIT_JSRUN_TAIL" CACHE INTERNAL "")
set(varvals ${varvals}  "|& grep -v \"IEEE_UNDERFLOW_FLAG\""  CACHE INTERNAL "")

#cprnc line
set(varnames ${varnames}  "SUMMIT_JSRUN_SERIAL"  CACHE INTERNAL "")
set(varvals ${varvals}  "jsrun -n 1 -r 1 -l cpu-cpu -b packed:1 -d plane:1 -a1 -c7 -g0" CACHE INTERNAL "")

#only for p9 runs
#set(varnames ${varnames}  "SUMMIT_OMP_OPTIONS" CACHE INTERNAL "")
#set(varvals ${varvals}  "export OMP_PLACES=threads\; export OMP_PROC_BIND=close\;export OMP_NUM_THREADS=4"  CACHE INTERNAL "")













