# Default guess at KOKKOS_OPTIONS. These can be overridden by
# appending -D${OPTION_NAME}=Off. The CMAKE_CXX_COMPILER will
# be handled automatically by CIME unless explicitly set in
# KOKKOS_OPTIONS. CMAKE_INSTALL_PREFIX will be handled automatically
# by CIME always.
set(KOKKOS_OPTIONS "-DKokkos_ENABLE_SERIAL=On")
if (compile_threaded)
  string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_OPENMP=On")
endif()

# Unless told otherwise, set has_contiguous to FALSE
set(HAS_F2008_CONTIGUOUS "FALSE")

# By default, link with CXX
set(E3SM_LINK_WITH_FORTRAN "FALSE")

# Do not use any automatic settings for these
set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_DEBUG "")
set(CMAKE_Fortran_FLAGS_RELEASE "")
set(CMAKE_C_FLAGS "")
set(CMAKE_C_FLAGS_DEBUG "")
set(CMAKE_C_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_Fortran_FORMAT_FIXED_FLAG "")
set(CMAKE_Fortran_FORMAT_FREE_FLAG "")


