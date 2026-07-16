# Water isotope physics process.
#
# This directory is only entered when EAMXX_ENABLE_WATER_ISOTOPES is ON, 
# which in turn guarantees EAMXX_ENABLE_WATER_TRACERS is ON
# (isotopes are layered on top of the water tracer infrastructure). Following
# the EAMxx one-library-per-process pattern, the water isotope process gets its
# own library.

add_library(water_isotopes
  eamxx_water_isotopes_process_interface.cpp
)

# Public compile definition so the registration guard in register_physics.hpp
# becomes active.
target_compile_definitions(water_isotopes PUBLIC EAMXX_HAS_WATER_ISOTOPES)

# Expose this directory
target_include_directories(water_isotopes PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(water_isotopes PUBLIC scream_share)
target_link_libraries(water_isotopes PRIVATE water_tracers)

# Link into the eamxx_physics INTERFACE target.
if (TARGET eamxx_physics)
  target_link_libraries(eamxx_physics INTERFACE water_isotopes)
endif()
