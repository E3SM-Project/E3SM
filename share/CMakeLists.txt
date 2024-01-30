#===============================================================================
# The CMake system for csm_share. This is mostly copied from the CMakeLists.txt
# for the main E3SM.
#===============================================================================

cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0057 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0079 NEW) # Remove once scorpio in a better state
set(CMAKE_CXX_STANDARD 17)

# We need to set the compilers *before* calling `project`.
# The only way to get the compiler name, is to load Macros.cmake
# However, we do *not* want to pollute the environment with other
# vars coming from Macros.cmake, so we encapsulate its inclusion
# in a new scope.
# Additionally, we also set CMAKE_BUILD_TYPE=DEBUG if Macros.cmake
# contains DEBUG set to true
function(set_compilers_e3sm)
  # Grab CXX compiler from CIME
  include(${CASEROOT}/Macros.cmake)

  if (DEBUG)
    set(CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "build type, inferred from ${DEBUG}")
  else()
    set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "build type, inferred from ${DEBUG}")
  endif()

  set(CMAKE_CXX_COMPILER     ${CMAKE_CXX_COMPILER}     CACHE STRING "The CXX compiler")
  set(CMAKE_C_COMPILER       ${CMAKE_C_COMPILER}       CACHE STRING "The C   compiler")
  set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} CACHE STRING "The Fortran compiler")
endfunction()

set_compilers_e3sm()

project(CSM_SHARE C CXX Fortran)

# Any changes to SourceMods will require us to reconfigure
file(GLOB COMPONENT_SOURCE_MOD_DIRS "${CASEROOT}/SourceMods/src.*")
foreach(COMPONENT_SOURCE_MOD_DIR IN LISTS COMPONENT_SOURCE_MOD_DIRS)
  set_property(
    DIRECTORY
    APPEND
    PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${COMPONENT_SOURCE_MOD_DIR})
endforeach()

# Include function definitions
include(${SRCROOT}/components/cmake/cmake_util.cmake)
include(${SRCROOT}/components/cmake/build_mpas_model.cmake)
include(${SRCROOT}/components/cmake/build_eamxx.cmake)
include(${SRCROOT}/components/cmake/build_model.cmake)

# Set up CMAKE_MODULE_PATH so any component can use E3SM
# and CIME cmake modules if they want.
list(APPEND CMAKE_MODULE_PATH ${SRCROOT}/components/cmake/modules)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/CIME/non_py/src/CMake)

if (COMP_INTERFACE STREQUAL "moab")
  set(USE_MOAB True)
  set(CPPDEFS "${CPPDEFS} -DHAVE_MOAB")
endif()

if (USE_MOAB)
  find_package(MOAB REQUIRED)
endif()

set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(BUILDCONF ${CASEROOT}/Buildconf)

# Set global targets
if (NOT TARGET genf90)
  add_custom_target(genf90
    DEPENDS ${CIMEROOT}/CIME/non_py/externals/genf90/genf90.pl)
endif()

# Build CSM_share
build_model("csm_share" "csm_share")
