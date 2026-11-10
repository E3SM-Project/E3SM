# CMake initial cache file for the ghci cuda image (gnu + CUDA)
#
# The GPU arch (Kokkos_ARCH_XYZ and CMAKE_CUDA_ARCHITECTURES) is NOT set here,
# since it depends on the GPU of the runner. Pass it via cmake options (see test_homme.py --cmake-opt).
include (${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

# CUDA compilation requires nvcc_wrapper (from Kokkos) as the C++ backend compiler of the MPI wrapper.
# NOTE: this cannot be done here, since env changes in this file only affect the cmake configure step,
# not the later build step. Before configuring/building, set (e.g., in the CI workflow):
#   NVCC_WRAPPER_DEFAULT_COMPILER=g++
#   OMPI_CXX=<path to externals/ekat/extern/kokkos/bin/nvcc_wrapper>   (if using OpenMPI)
#   MPICH_CXX=<path to externals/ekat/extern/kokkos/bin/nvcc_wrapper>  (if using MPICH)

SET (HOMMEXX_EXEC_SPACE CUDA CACHE STRING "")

# Our runners' MPI is not cuda-aware, so stage comms buffers through host memory
# TODO: set to ON (or remove) once we have cuda-aware MPI
SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")

SET (Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
# Compose (always) needs the Kokkos Serial backend. Kokkos disables it by default if another
# host backend (e.g., OpenMP) is on, so we must enable it explicitly.
SET (Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
SET (Kokkos_ENABLE_CUDA ON CACHE BOOL "")
SET (Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
SET (Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

SET (CXXLIB_SUPPORTED_CACHE FALSE CACHE BOOL "")

SET (ENABLE_OPENMP OFF CACHE BOOL "")
SET (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
SET (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
SET (BLAS_LIBRARIES   "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "" FORCE)
SET (LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "" FORCE)
