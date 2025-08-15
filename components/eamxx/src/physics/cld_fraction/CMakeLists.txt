set(CLDFRAC_SRCS
  eamxx_cld_fraction_process_interface.cpp
  cld_fraction.cpp
)

set(CLDFRAC_HEADERS
  eamxx_cld_fraction_process_interface.hpp
  cld_fraction_functions.hpp
  cld_fraction_main_impl.hpp
)

add_library(cld_fraction ${CLDFRAC_SRCS})
target_compile_definitions(cld_fraction PUBLIC EAMXX_HAS_CLD_FRACTION)
target_link_libraries(cld_fraction PUBLIC physics_share scream_share)

option (EAMXX_ENABLE_CLDFRAC_LAPIS_SUPPORT "Whether CldFraction should support a LAPIS-generated ML emulator" OFF)
if (EAMXX_ENABLE_CLDFRAC_LAPIS_SUPPORT)
  target_compile_definitions(cld_fraction PRIVATE EAMXX_ENABLE_CLDFRAC_LAPIS_SUPPORT)

  # Build LAPIS-generated cpp file
  add_library(cldfrac_ml_lapis cldfrac_ml_lapis.cpp)
  target_link_libraries (cldfrac_ml_lapis PUBLIC Kokkos::kokkos)

  target_link_libraries(cld_fraction PRIVATE cldfrac_ml_lapis)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE cld_fraction)
endif()

# Cloud fraction does not yet have a set of unit tests or a BFB test comparing with the F90
# code.
# The cloud fraction stand alone test, in the /tests/ directory covers a range of property
# tests.
# The BFB test will require some amount of work, as the specific ice cloud fraction scheme
# here is just an option in a more involved F90 subroutine.  Furture work would be to
# include a BFB test with the F90 code, located in /components/eam/src/physics/cam/cldfrc2m.F90
