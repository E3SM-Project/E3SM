set(CLDFRAC_SRCS
  atmosphere_cld_fraction.cpp
  cld_fraction.cpp 
)

set(CLDFRAC_HEADERS
  atmosphere_cld_fraction.hpp
  cld_fraction_functions.hpp
  cld_fraction_main_impl.hpp
)

add_library(cld_fraction ${CLDFRAC_SRCS})
target_include_directories(cld_fraction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(cld_fraction physics_share scream_share)
target_compile_options(cld_fraction PUBLIC)

# Cloud fraction does not yet have a set of unit tests or a BFB test comparing with the F90
# code.
# The cloud fraction stand alone test, in the /tests/ directory covers a range of property
# tests.
# The BFB test will require some amount of work, as the specific ice cloud fraction scheme
# here is just an option in a more involved F90 subroutine.  Furture work would be to
# include a BFB test with the F90 code, located in /components/eam/src/physics/cam/cldfrc2m.F90
