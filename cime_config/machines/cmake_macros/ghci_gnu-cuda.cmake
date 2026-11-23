# Pick up the generic gnu settings (CIME would not load them automatically for this compiler name)
include (${CMAKE_CURRENT_LIST_DIR}/gnu.cmake)

string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fcheck=bounds")

string (APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON")

# Change to ON if/when our runners feature cuda-aware MPI
set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "")

# Cap the Hommexx team size to 8 warps. With the default (16 warps), some homme kernels
# need too many registers, and abort at runtime with "Requested too large team size.
# Requested: 16, Maximum: 8". Same cap as in ghci-hip.cmake.
SET (HOMMEXX_CUDA_MAX_WARP_PER_TEAM 8 CACHE STRING "")

