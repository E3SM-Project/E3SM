# CMake initial cache file for the ghci gnu image
include (${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

SET (Kokkos_ARCH_ZEN2 ON CACHE BOOL "")

# Kokkos CPU parallelism: threads via OpenMP
SET (ENABLE_OPENMP OFF CACHE BOOL "")
SET (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
SET (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")
SET (Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
# Compose (always) needs the Kokkos Serial backend. Kokkos disables it by default if another
# host backend (e.g., OpenMP) is on, so we must enable it explicitly.
SET (Kokkos_ENABLE_SERIAL ON CACHE BOOL "")

SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
SET (BLAS_LIBRARIES   "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "" FORCE)
SET (LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "" FORCE)
