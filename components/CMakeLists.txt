#===============================================================================
#
# Common CMakeLists.txt: a framework for building all CIME components and more
#
# This is a port of cime/CIME/Tools/Makefile. As more components are ported to
# CMake, the directory level of this file will rise to the top-level directory.
#
# We will prefer space-separated strings over lists
#
#===============================================================================

cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0057 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0079 NEW) # Remove once scorpio in a better state
set(CMAKE_CXX_STANDARD 17)

# Turn on wrapper
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "python3 ${CASEROOT}/Tools/e3sm_compile_wrap.py ")
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK    "python3 ${CASEROOT}/Tools/e3sm_compile_wrap.py ")

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
    set(E3SM_DEFAULT_BUILD_TYPE "DEBUG" CACHE STRING "Default build type, inferred from ${DEBUG}")
  else()
    set(E3SM_DEFAULT_BUILD_TYPE "RELEASE" CACHE STRING "Default build type, inferred from ${DEBUG}")
  endif()

  set(CMAKE_CXX_COMPILER     ${CMAKE_CXX_COMPILER}     CACHE STRING "The CXX compiler")
  set(CMAKE_C_COMPILER       ${CMAKE_C_COMPILER}       CACHE STRING "The C   compiler")
  set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} CACHE STRING "The Fortran compiler")

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

# If no CMAKE_BUILD_TYPE/CMAKE_CONFIGURATION_TYPES are provided,
# then set the build type according to the DEBUG variable found in Macros.cmake

# Why do we CMAKE_BUILD_TYPE to be set? Glad you asked. Read on.
# Some subfolders might need it (or even try to set it based on
# Macros.cmake anyways). However, not setting it GLOBALLY
# has the risk that cmake may treat targets differently in different subfolders.
# One BIG exapmle: the yaml-cpp target appends a 'd' to its file name when
# the build type is DEBUG (i.e., the lib name is libyaml-cppd.a).
# If yaml-cpp is built in a subdir that locally sets CMAKE_BUILD_TYPE=DEBUG,
# cmake will create the yaml-cpp target to point to libyaml-cppd.a.
# However, in other directories where CMAKE_BUILD_TYPE is not
# set (or not set to DEBUG), cmake will try to  link libyaml-cpp.a
# Bottom line: just like CMAKE_<LANG>_COMPILER, you SHOULD set
# CMAKE_BUILD_TYPE at most once, and do it BEFORE calling `project`.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)

  message(STATUS "Setting build type to '${E3SM_DEFAULT_BUILD_TYPE}', as none was specified. This default comes from what CIME set in the DEBUG variable")
  set(CMAKE_BUILD_TYPE "${E3SM_DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)

  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "DEBUG" "RELEASE")
endif()

project(E3SM C CXX Fortran)

if (USE_CUDA)
  set(CMAKE_CUDA_COMPILER_FORCED True)
  enable_language(CUDA)
elseif (USE_HIP)
  set(CMAKE_HIP_COMPILER_FORCED True)
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
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_util.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_mpas_model.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_eamxx.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_model.cmake)

# Set up CMAKE_MODULE_PATH so any component can use E3SM
# and CIME cmake modules if they want.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/CIME/non_py/src/CMake)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Find dependencies
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/find_dep_packages.cmake)

# Scream manages its own flags
build_eamxx()

set(BUILDCONF ${CASEROOT}/Buildconf)

build_mpas_models()

# Set global cmake settings
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# Set global targets
if (NOT TARGET genf90)
  add_custom_target(genf90
    DEPENDS ${CIMEROOT}/CIME/non_py/externals/genf90/genf90.pl)
endif()

# Build E3SM components
set(IDX 0)
set(COMP_CLASSES cpl atm lnd ice ocn rof glc wav iac esp)
set(SKIP_COMPS "scream" "mpaso" "mpassi" "mali")
foreach(COMP_NAME IN LISTS COMP_NAMES)
  list(GET COMP_CLASSES ${IDX} COMP_CLASS)
  if (NOT COMP_CLASS STREQUAL "cpl" AND NOT COMP_NAME IN_LIST SKIP_COMPS)
    message("Found component ${COMP_CLASS} model '${COMP_NAME}'")
    add_subdirectory(cmake/${COMP_CLASS})
  endif()

  math(EXPR IDX "${IDX} + 1")
endforeach()

# Build the E3SM coupler and exe
add_subdirectory(cmake/cpl)
