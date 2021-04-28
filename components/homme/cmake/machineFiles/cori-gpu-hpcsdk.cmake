# CMake initial cache file for Cori gpu nodes (named cori-gpu -- small testbed with V100's and SKX cpus)

SET(HOMMEXX_EXEC_SPACE CUDA CACHE STRING "")
#SET(HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")
#SET(HOMMEXX_CUDA_MAX_WARP_PER_TEAM "16" CACHE STRING  "")

SET(NETCDF_DIR $ENV{NETCDF_DIR} CACHE FILEPATH "")
SET(PNETCDF_DIR $ENV{PARALLEL_NETCDF_DIR} CACHE FILEPATH "")
SET(HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")

SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

SET(ENABLE_CUDA FALSE CACHE BOOL "")

#SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")

#SET(HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")

SET(USE_TRILINOS OFF CACHE BOOL "")

SET(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA ON CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
SET(Kokkos_ARCH_VOLTA70 ON CACHE BOOL "")
#SET(Kokkos_ARCH_SKX ON CACHE BOOL "") # ndk -- not sure why it didn't like this on as well pgc++-Error-Unknown switch: -mtune=skylake-avx512
#SET(Kokkos_ENABLE_CUDA_UVM ON CACHE BOOL "") #ndk -- get build error when I tried UVM
SET(Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA_ARCH_LINKING OFF CACHE BOOL "") #ndk -- need this to get around link error with fortran (-arch=sm_70) 

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "")
# ndk -- need to set OMPI_CXX env variable, but can't seem to set it here
#SET(ENV{OMPI_CXX} "/global/cscratch1/sd/ndk/wacmy/m54-apr16/components/homme/../../externals/kokkos/bin/nvcc_wrapper" CACHE FILEPATH "")
#SET(ENV{OMPI_CXX} "/global/cscratch1/sd/ndk/wacmy/m54-apr16/components/homme/../../externals/kokkos/bin/nvcc_wrapper")
#message("ndk ENV{OMPI_CXX}=$ENV{OMPI_CXX}")
# ndk -- also need to change default host_compiler in externals/kokkos/bin/nvcc_wrapper which can be done via:
#setenv NVCC_WRAPPER_DEFAULT_COMPILER "pgc++"


#SET(E3SM_KOKKOS_PATH "$KOKKOS_DIR" CACHE BOOL "")
#SET(E3SM_KOKKOS_PATH "/usr/common/software/sles15_cgpu/kokkos/3.2.00" CACHE BOOL "")

SET(ENABLE_OPENMP OFF CACHE BOOL "")
SET(ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
SET(ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

#SET(HOMME_TESTING_PROFILE "dev" CACHE STRING "")

SET(USE_NUM_PROCS 4 CACHE STRING "")
