string (APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_SPR=ON")

# Help find_package to find the correct BLAS (i.e., the MKL installations)
set (BLA_VENDOR "Intel10_64lp" CACHE STRING "")
