string (APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON")

# Change to ON if/when our runners feature cuda-aware MPI
set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "")
