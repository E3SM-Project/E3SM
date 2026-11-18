include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(EKAT_MPI_EXTRA_ARGS "${EKAT_MPI_EXTRA_ARGS} --gpus-per-task=1" CACHE STRING "" FORCE)

set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)
set(HOMME_QUAD_PREC FALSE CACHE BOOL "") # nvidia does not support QUAD (real*16)

option(SCREAM_SMALL_KERNELS "Use small, non-monolothic kokkos kernels for ALL components that support them" ON)
