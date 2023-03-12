include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

set(BLAS_LIBRARIES /opt/conda/lib/libblas.so CACHE STRING "")
set(LAPACK_LIBRARIES /opt/conda/lib/liblapack.so CACHE STRING "")
set(SCREAM_INPUT_ROOT "/storage/inputdata/" CACHE STRING "")