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

  if (MPILIB STREQUAL "mpi-serial")
    set(CC ${SCC})
    set(FC ${SFC})
    set(CXX ${SCXX})
  else()
    set(CC ${MPICC})
    set(FC ${MPIFC})
    set(CXX ${MPICXX})
  endif()

  set(CMAKE_CXX_COMPILER     ${CXX} CACHE STRING "The CXX compiler")
  set(CMAKE_C_COMPILER       ${CC}  CACHE STRING "The C   compiler")
  set(CMAKE_Fortran_COMPILER ${FC}  CACHE STRING "The Fortran compiler")

  # USE_CUDA or USE_HIP is set through Macros.cmake
  # For instance: cime_config/machines/cmake_macros/gnugpu_summit.cmake
  # If it exists, then set parent's scope to true; otherwise to false
  # At this point, we use either CUDA or HIP.
  # Revisit as needed for future systems.
  if (USE_CUDA)
    set(USE_CUDA TRUE PARENT_SCOPE)
  elseif (USE_HIP)
    set(USE_HIP TRUE PARENT_SCOPE)
  else()
    set(USE_CUDA FALSE PARENT_SCOPE)
    set(USE_HIP FALSE PARENT_SCOPE)
  endif()
endfunction()

set_compilers_e3sm()

project(CSM_SHARE C CXX Fortran)

if(USE_CUDA)
  enable_language(CUDA)
elseif(USE_HIP)
  enable_language(HIP)
endif()

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

set(CMAKE_VERBOSE_MAKEFILE TRUE)

# We do want CMAKE_BUILD_TYPE to be set, but we do NOT want CMake to
# decide what optimization flags to append, based on build type,
# for components who rely on CIME for build flags, so make all the following empty.
# JGF: I think we can remove this once proper CMake names are being used in the
# macros.
set (CMAKE_C_FLAGS_RELEASE "")
set (CMAKE_CXX_FLAGS_RELEASE "")
set (CMAKE_Fortran_FLAGS_RELEASE "")

set (CMAKE_C_FLAGS_DEBUG "")
set (CMAKE_CXX_FLAGS_DEBUG "")
set (CMAKE_Fortran_FLAGS_DEBUG "")

set(BUILDCONF ${CASEROOT}/Buildconf)

# Set global targets
if (NOT TARGET genf90)
  add_custom_target(genf90
    DEPENDS ${CIMEROOT}/CIME/non_py/externals/genf90/genf90.pl)
endif()

# Build CSM_share
build_model("csm_share" "csm_share")
