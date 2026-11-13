# CMake initial cache file for the ghci hip image (gnu + HIP)
#
# The GPU arch (Kokkos_ARCH_AMD_GFXXYZ) is NOT set here, since it depends on the GPU of the runner.
# Pass it via cmake options (see test_homme.py --cmake-opt).
include (${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

# Only build (and test) the kokkos dycores: the F90 executables run on the host, so their tests
# are already covered by the CPU CI jobs, and they would take most of the test time on the GPU runners.
# NOTE: with HOMMEXX_BFB_TESTING=ON, the F90 theta/preqx dycores are turned back on (for the cxx_vs_f90 tests).
SET (BUILD_HOMME_THETA OFF CACHE BOOL "")
SET (BUILD_HOMME_PREQX OFF CACHE BOOL "")
SET (BUILD_HOMME_SWEQX OFF CACHE BOOL "")
SET (BUILD_HOMME_TOOL  OFF CACHE BOOL "")

# HIP compilation requires a HIP compiler (hipcc) as the C++ backend compiler of the MPI wrapper.
# NOTE: this cannot be done here, since env changes in this file only affect the cmake configure step,
# not the later build step. Before configuring/building, set (e.g., in the CI workflow):
#   OMPI_CXX=hipcc   (if using OpenMPI)
#   MPICH_CXX=hipcc  (if using MPICH)

SET (HOMMEXX_EXEC_SPACE HIP CACHE STRING "")

# Our runners' MPI is not GPU-aware, so stage comms buffers through host memory
# TODO: set to ON (or remove) once we have GPU-aware MPI
SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")

SET (Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
# Compose (always) needs the Kokkos Serial backend. Kokkos disables it by default if another
# host backend (e.g., OpenMP) is on, so we must enable it explicitly.
SET (Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
SET (Kokkos_ENABLE_HIP ON CACHE BOOL "")
SET (Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

SET (CXXLIB_SUPPORTED_CACHE FALSE CACHE BOOL "")

SET (ENABLE_OPENMP OFF CACHE BOOL "")
SET (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
SET (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
SET (BLAS_LIBRARIES   "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "" FORCE)
SET (LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "" FORCE)
