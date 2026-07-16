# Auxiliary tracers for EAMxx.
#
# This folder groups the "auxiliary tracer" physics processes. Currently,
# there are two options: passive water tracers and, built on top of them, 
# water isotopes. Additional tracer classes could be added in an analagous
# way. The build is governed by these CMake options:
#
#   EAMXX_ENABLE_WATER_TRACERS  - build the passive water tracer process
#   EAMXX_ENABLE_WATER_ISOTOPES - build the water isotope process
#
# Water isotopes are implemented on top of the water tracer infrastructure, so
# they can ONLY be built when water tracers are also enabled. If isotopes are
# requested without tracers we auto-enable tracers and warn, rather than build
# an inconsistent configuration.

option(EAMXX_ENABLE_WATER_TRACERS  "Whether to build the EAMxx water tracer process"  OFF)
option(EAMXX_ENABLE_WATER_ISOTOPES "Whether to build the EAMxx water isotope process" OFF)

# Enforce the dependency: isotopes require tracers.
if (EAMXX_ENABLE_WATER_ISOTOPES AND NOT EAMXX_ENABLE_WATER_TRACERS)
  message(STATUS "WARNING: EAMXX_ENABLE_WATER_ISOTOPES=ON requires EAMXX_ENABLE_WATER_TRACERS; "
                 "auto-enabling water tracers.")
  set(EAMXX_ENABLE_WATER_TRACERS ON CACHE BOOL "Whether to build the EAMxx water tracer process" FORCE)
endif()

# Each process lives in its own sibling subdirectory. Descend
# into each only when the corresponding process is enabled. water_isotopes must
# be added after water_tracers, since it depends on the water_tracers target.
if (EAMXX_ENABLE_WATER_TRACERS)
  add_subdirectory(water_tracers)
endif()
if (EAMXX_ENABLE_WATER_ISOTOPES)
  add_subdirectory(water_isotopes)
endif()
