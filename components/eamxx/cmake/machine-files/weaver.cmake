# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/weaver.cmake)
set (BLAS_LIBRARIES /ascldap/users/projects/e3sm/scream/libs/openblas/install/weaver/gcc/8.5.0/lib/libopenblas.so CACHE STRING "")
set (LAPACK_LIBRARIES /ascldap/users/projects/e3sm/scream/libs/openblas/install/weaver/gcc/8.5.0/lib/libopenblas.so CACHE STRING "")
set(SCREAM_INPUT_ROOT "/home/projects/e3sm/scream/data" CACHE STRING "")
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)
set(HOMMEXX_CUDA_MAX_WARP_PER_TEAM 8 CACHE STRING "")
