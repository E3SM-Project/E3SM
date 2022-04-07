#===============================================================================
#
# Common CMakeLists.txt: a framework for building all CIME components and more
#
# This is a port of cime/scripts/Tools/Makefile. As more components are ported to
# CMake, the directory level of this file will rise to the top-level directory.
#
# We will prefer space-separated strings over lists
#
#===============================================================================

# bmpersch
# - Read modern cmake docs, use modern features
# - Use find_package for trilinos and other TPLS

cmake_minimum_required(VERSION 3.9)
cmake_policy(SET CMP0057 NEW)

# Store caseroot in the cache, so that, if cmake re-runs,
# we can still access ${CASEROOT}/Macros.cmake
set (CASEROOT "" CACHE STRING "The case root directory.")

# Turn on wrapper
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "python3 ${CASEROOT}/Tools/e3sm_compile_wrap.py ")
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK    "python3 ${CASEROOT}/Tools/e3sm_compile_wrap.py ")

# We need to set the compilers *before* calling `project`.
# The only way to get the compiler name, is to load Macros.cmake
# However, we do *not* want to pollut the environment with other
# vars coming from Macros.cmake, so we encapsulate its inclusion
# in a new scope.
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
  # USE_CUDA is set through Macros.cmake / config_compilers.xml
  # If it exists, then set parent's scope to true; otherwise to false
  if (USE_CUDA)
    set(USE_CUDA TRUE PARENT_SCOPE)
  else()
    set(USE_CUDA FALSE PARENT_SCOPE)
  endif()
endfunction()
set_compilers_e3sm()

project(E3SM C CXX Fortran)

if(USE_CUDA)
  enable_language(CUDA)
endif()

# Include function definitions
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_util.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_mpas_model.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_model.cmake)

set(BUILDCONF ${CASEROOT}/Buildconf)

# Do any MPAS builds first since they manage their own flags
build_mpas_models()

# Set global cmake settings
set(CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# Set global targets
add_custom_target(genf90
  DEPENDS ${CIMEROOT}/CIME/non_py/externals/genf90/genf90.pl)

# Build E3SM components
set(IDX 0)
set(COMP_CLASSES cpl atm lnd ice ocn rof glc wav iac esp)
set(MPAS_COMPS "mpaso" "mpassi" "mali")
foreach(COMP_NAME IN LISTS COMP_NAMES)
  list(GET COMP_CLASSES ${IDX} COMP_CLASS)
  if (NOT COMP_CLASS STREQUAL "cpl" AND NOT COMP_NAME IN_LIST MPAS_COMPS)
    message("Found component ${COMP_CLASS} model '${COMP_NAME}'")
    add_subdirectory(cmake/${COMP_CLASS})
  endif()

  math(EXPR IDX "${IDX} + 1")
endforeach()

add_subdirectory(cmake/cpl)
