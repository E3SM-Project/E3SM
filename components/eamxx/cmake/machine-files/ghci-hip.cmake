# Common settings for our ghci images
include(${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

# Set Fortran flags
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch" CACHE STRING "Fortran compiler flags" FORCE)

# Set the path to BLAS/LAPACK libraries
set(BLAS_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "Path to BLAS library" FORCE)
set(LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "Path to LAPACK library" FORCE)

# Set SCREAM_MACHINE
set(SCREAM_MACHINE ghci-hip CACHE STRING "")

# Enable HIP in kokkos
# NOTE: the GPU arch (Kokkos_ARCH_AMD_GFXXXX) is NOT set here, since it depends on the GPU of
#       the runner. It is passed via cmake options (see the get-hip-arch action).
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)

# Same as on frontier
set(CMAKE_CXX_FLAGS "-Wno-mismatched-tags -munsafe-fp-atomics -fno-gpu-rdc" CACHE STRING "" FORCE)

# Our runners' MPI is not gpu-aware
option(SCREAM_MPI_ON_DEVICE "Whether to use device pointers for MPI calls" OFF)

# We have 2 GPUs/node, and we run a SINGLE build per node, so we can fit 2 ranks there
set(SCREAM_TEST_MAX_RANKS 2 CACHE STRING "Upper limit on ranks for mpi tests")
