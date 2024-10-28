#module restore
#module load oneapi/eng-compiler/2022.12.30.005
#module load intel_compute_runtime/release/agama-devel-627
#module load spack cmake
#module list



SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")
SET(HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

#temp hack
SET(HOMME_USE_KOKKOS TRUE CACHE BOOL "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")

#set(KOKKOS_HOME "/home/onguba/kokkos-build/jan03-2024/install" CACHE STRING "")
#set(E3SM_KOKKOS_PATH ${KOKKOS_HOME} CACHE STRING "")

SET(USE_TRILINOS OFF CACHE BOOL "")

SET(SYCL_BUILD TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")

SET(CMAKE_CXX_STANDARD 17)

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "")

# -fsycl-link-huge-device-code for theta to get build
SET(SYCL_COMPILE_FLAGS "-std=c++17 -fsycl -fsycl-device-code-split=per_kernel -fno-sycl-id-queries-fit-in-int -fsycl-unnamed-lambda")
SET(SYCL_LINK_FLAGS "-fsycl -fsycl-link-huge-device-code -fsycl-device-code-split=per_kernel -fsycl-targets=spir64")

SET(ADD_Fortran_FLAGS "-fc=ifx -O3 -DNDEBUG -DCPRINTEL -g" CACHE STRING "")
SET(ADD_C_FLAGS "-O3 -DNDEBUG " CACHE STRING "")

SET(ADD_CXX_FLAGS "-std=c++17 -O3 -DNDEBUG ${SYCL_COMPILE_FLAGS}" CACHE STRING "")
SET(ADD_LINKER_FLAGS "-O3 -DNDEBUG ${SYCL_LINK_FLAGS} -fortlib" CACHE STRING "")

set (ENABLE_OPENMP OFF CACHE BOOL "")
set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

set (HOMME_TESTING_PROFILE "dev" CACHE STRING "")

set (USE_NUM_PROCS 4 CACHE STRING "")

SET (USE_MPI_OPTIONS "--bind-to core" CACHE FILEPATH "")


