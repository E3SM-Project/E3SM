# CMake initial cache file for the ghci intel image
include (${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

SET (Kokkos_ARCH_SPR ON CACHE BOOL "")

SET (ENABLE_OPENMP OFF CACHE BOOL "")
SET (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
SET (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")
SET (Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
# Compose (always) needs the Kokkos Serial backend. Kokkos disables it by default if another
# host backend (e.g., OpenMP) is on, so we must enable it explicitly.
SET (Kokkos_ENABLE_SERIAL ON CACHE BOOL "")

# Use MKL for blas/lapack
SET (HOMME_USE_MKL TRUE CACHE BOOL "")

# Assuming this is the (clang-based) icpx compiler, this flag turns off annoying compiler warning
SET (ADD_CXX_FLAGS "-Wno-overriding-option" CACHE STRING "")
