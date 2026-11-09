# Common settings for our ghci images
include(${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

# Set Fortran flags
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch" CACHE STRING "Fortran compiler flags" FORCE)

# Set the path to BLAS/LAPACK libraries
set(BLAS_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "Path to BLAS library" FORCE)
set(LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "Path to LAPACK library" FORCE)

# Set SCREAM_MACHINE
set(SCREAM_MACHINE ghci-hip CACHE STRING "")

# Enable HIP in kokkos (the arch is passed in by the CI workflow, via Kokkos_ARCH_AMD_GFX* cmake config)
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)

# Our runners do not (yet) feature GPU-aware MPI
option(SCREAM_MPI_ON_DEVICE "Whether to use device pointers for MPI calls" OFF)

# Python tests rely on the venv in ghci.cmake, which is available on the hip image too
option (EAMXX_ENABLE_PYTHON "Whether to enable python interface from eamxx" ON)

# On gfx90a, a workgroup of 1024 threads leaves at most 128 registers per thread (VGPRs + AGPRs). With -DNDEBUG,
# ROCm 6.4.2 compiles some Hommexx kernels (e.g., DirkFunctorImpl::run_newton) with 132 registers, so launching
# them with 1024 threads (16 warps/team, the default) aborts with HSA_STATUS_ERROR_INVALID_ISA.
# Cap the team size to 8 warps (512 threads).
set (HOMMEXX_CUDA_MAX_WARP_PER_TEAM 8 CACHE STRING "Maximum number of warps per team on GPU")
