# Auxiliary tracers for EAMxx.
#
# This folder groups the "auxiliary tracer" physics processes. Currently
# there is one process, water_tracers, which accepts generic and isotope
# tracer entries in a single water_tracers_file. The build is governed by these
# CMake options:
#
#   EAMXX_ENABLE_WATER_TRACERS  - build the water tracer process
#   EAMXX_ENABLE_WATER_ISOTOPES - allow `type: isotope` entries in the
#                                 tracers file (isotope math is header-only
#                                 and always compiled into water_tracers;
#                                 this flag only gates acceptance and the
#                                 isotope-specific unit tests)
#
# Isotopes are implemented on top of the water tracer infrastructure, so they
# can ONLY be enabled when water tracers are also enabled. If isotopes are
# requested without tracers we auto-enable tracers and warn, rather than build
# an inconsistent configuration.

option(EAMXX_ENABLE_WATER_TRACERS  "Whether to build the EAMxx water tracer process"  OFF)
option(EAMXX_ENABLE_WATER_ISOTOPES "Whether to accept isotope entries in the EAMxx water tracer process" OFF)

# Enforce the dependency: isotopes require tracers.
if (EAMXX_ENABLE_WATER_ISOTOPES AND NOT EAMXX_ENABLE_WATER_TRACERS)
  message(STATUS "WARNING: EAMXX_ENABLE_WATER_ISOTOPES=ON requires EAMXX_ENABLE_WATER_TRACERS; "
                 "auto-enabling water tracers.")
  set(EAMXX_ENABLE_WATER_TRACERS ON CACHE BOOL "Whether to build the EAMxx water tracer process" FORCE)
endif()

if (EAMXX_ENABLE_WATER_TRACERS)
  add_subdirectory(water_tracers)
endif()
