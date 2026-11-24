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

# Cap the Hommexx team size to 8 warps (512 threads on AMD GPUs). With the default (16 warps, i.e., 1024 threads),
# optimized builds (-DNDEBUG) abort at runtime with HSA_STATUS_ERROR_INVALID_ISA on our MI210 (gfx90a) runners:
# with ROCm 6.4.2, some Hommexx kernels (e.g., Homme::DirkFunctorImpl::run_newton) end up needing 132 registers per
# thread, more than the 128 available when launching 1024-thread workgroups on gfx90a.
# NOTE: Frontier (ROCm 6.2.4/7.2.0) uses 16, so this may no longer be needed with another ROCm version in the image.
set (HOMMEXX_CUDA_MAX_WARP_PER_TEAM 8 CACHE STRING "Maximum number of warps per team on GPU")

