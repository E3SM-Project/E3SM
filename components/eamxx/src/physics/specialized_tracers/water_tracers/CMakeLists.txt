# Water tracer physics process.
#
# This directory is only entered when EAMXX_ENABLE_WATER_TRACERS is ON. 
# Following the EAMxx one-library-per-process pattern, the
# water tracer process gets its own library. The sibling water_isotopes library
# depends on this one.

add_library(water_tracers
  eamxx_water_tracers_process_interface.cpp
)

# Public compile definition so the registration guard in register_physics.hpp
# becomes active.
target_compile_definitions(water_tracers PUBLIC EAMXX_HAS_WATER_TRACERS)

# Expose this directory to allow includes of the
# process interface headers.
target_include_directories(water_tracers PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link required EAMxx libraries.
target_link_libraries(water_tracers scream_share)

# Link into the eamxx_physics INTERFACE target.
if (TARGET eamxx_physics)
  target_link_libraries(eamxx_physics INTERFACE water_tracers)
endif()
