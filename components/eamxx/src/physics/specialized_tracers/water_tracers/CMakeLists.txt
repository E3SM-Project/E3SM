# Water tracer physics process.
#
# This directory is only entered when EAMXX_ENABLE_WATER_TRACERS is ON. There
# is a single process library: isotope-specific state and math (the sibling
# water_isotopes/ directory) are always compiled in -- they are header-only,
# pure Kokkos-callable math, so there is no cost to including them even when
# EAMXX_HAS_WATER_ISOTOPES is not defined. The macro only gates whether the
# process accepts `type: isotope` entries in the tracers file (see
# eamxx_water_tracers_process_interface.cpp).

add_library(water_tracers
  eamxx_water_tracers_process_interface.cpp
  eamxx_water_tracers_registry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../water_isotopes/eamxx_water_isotopes.cpp
)

# Public compile definition so the registration guard in register_physics.hpp
# becomes active.
target_compile_definitions(water_tracers PUBLIC EAMXX_HAS_WATER_TRACERS)

# EAMXX_ENABLE_WATER_ISOTOPES only gates isotope-entry acceptance in the
# registry/process; the isotope math itself is always compiled in.
if (EAMXX_ENABLE_WATER_ISOTOPES)
  target_compile_definitions(water_tracers PUBLIC EAMXX_HAS_WATER_ISOTOPES)
endif()

# Expose this directory and the sibling water_isotopes/ directory, so the
# process interface header can include eamxx_water_isotopes.hpp.
target_include_directories(water_tracers PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../water_isotopes
)

# Link required EAMxx libraries. The registry parses the tracers file with
# ekat's YAML parser; the process uses the common physics functions.
target_link_libraries(water_tracers PUBLIC scream_share eamxx_physics_share)
target_link_libraries(water_tracers PRIVATE ekat::YamlParser)

# Link into the eamxx_physics INTERFACE target.
if (TARGET eamxx_physics)
  target_link_libraries(eamxx_physics INTERFACE water_tracers)
endif()

# Unit tests (skipped in library-only builds).
if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
