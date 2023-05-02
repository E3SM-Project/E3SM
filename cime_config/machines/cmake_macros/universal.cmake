set(SUPPORTS_CXX "FALSE")

# Default guess at KOKKOS_OPTIONS. These can be overridden by
# appending -D${OPTION_NAME}=Off. The CMAKE_CXX_COMPILER will
# be handled automatically by CIME unless explicitly set in
# KOKKOS_OPTIONS. CMAKE_INSTALL_PREFIX will be handled automatically
# by CIME always.
set(KOKKOS_OPTIONS "-DKokkos_ENABLE_SERIAL=On")
if (compile_threaded)
  string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_OPENMP=On")
endif()
