 #
# This is the interface between E3SM's new CMake-based build system and MPAS.
#
# The following CMake variables are expected to be defined:
#  * CORES : A list of CORES to build, comma-separated (e.g. "ocean,seaice,landice")
#  * Whatever CIME settings are setting to correctly resolve the ${CASEROOT}/Macros.cmake file
#    - COMPILER, DEBUG, MPILIB, MACH, OS
#

# Source CIME-generated Macros
include(${CASEROOT}/Macros.cmake)
# Load machine/compiler specific settings
set(COMPILER_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${COMPILER}.cmake)
set(MACHINE_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${MACH}.cmake)
set(PLATFORM_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${MACH}.${COMPILER}.cmake)
set(TRY_TO_LOAD ${COMPILER_SPECIFIC_DEPENDS} ${MACHINE_SPECIFIC_DEPENDS} ${PLATFORM_SPECIFIC_DEPENDS})
foreach(ITEM IN LISTS TRY_TO_LOAD)
  if (EXISTS ${ITEM})
    include(${ITEM})
  endif()
endforeach()

#
# General setup
#

if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()

# Set up CPPDEFS
set(FILE_OFFSET "-DOFFSET64BIT")
if (CPPDEFS)
  separate_arguments(CPPDEFS UNIX_COMMAND "${CPPDEFS}")
endif()
list(APPEND CPPDEFS "-DMPAS_NO_LOG_REDIRECT" "-DMPAS_NO_ESMF_INIT" "-DMPAS_ESM_SHR_CONST" "-DMPAS_PERF_MOD_TIMERS" "${MODEL_FORMULATION}" "${FILE_OFFSET}" "${ZOLTAN_DEFINE}" "-D_MPI" "-DMPAS_NAMELIST_SUFFIX=${NAMELIST_SUFFIX}" "-DMPAS_EXE_NAME=${EXE_NAME}")
if (PIO_VERSION STREQUAL "2")
  list(APPEND CPPDEFS "-DUSE_PIO2")
  if (TARGET adios2::adios2)
    list(APPEND CPPDEFS "-DHAVE_ADIOS")
  endif()
endif()

if (DEBUG)
  list(APPEND CPPDEFS "-DMPAS_DEBUG")
endif()
if (compile_threaded)
  list(APPEND CPPDEFS "-DMPAS_OPENMP")
endif()

set(INCLUDES "${INSTALL_SHAREDPATH}/include" "${CMAKE_CURRENT_SOURCE_DIR}/external/ezxml" "${CMAKE_BINARY_DIR}/framework" "${CMAKE_BINARY_DIR}/operators")

if (USE_KOKKOS)
  include(${INSTALL_SHAREDPATH}/kokkos_generated_settings.cmake)
  string (REPLACE ";" " " KOKKOS_CXXFLAGS_STR "${KOKKOS_CXXFLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${KOKKOS_CXXFLAGS_STR}")
endif()

# Include custom cmake libraries used for mpas
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_utils.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_core.cmake)

# Add tools
add_subdirectory(tools)

# Gather sources that are needed for all cores into "common" library

set(COMMON_RAW_SOURCES external/ezxml/ezxml.c)
include(${CMAKE_CURRENT_SOURCE_DIR}/framework/framework.cmake)
if (COMP_INTERFACE STREQUAL "moab")
  list(APPEND COMMON_RAW_SOURCES framework/mpas_moabmesh.F)
  list(APPEND INCLUDES ${MOAB_INCLUDE_DIRS})
  list(APPEND CPPDEFS "-DHAVE_MOAB")
endif()
include(${CMAKE_CURRENT_SOURCE_DIR}/operators/operators.cmake)

add_library(common)
target_compile_definitions(common PRIVATE ${CPPDEFS})
target_include_directories(common PRIVATE ${INCLUDES})
target_link_libraries(common PUBLIC spio csm_share)

genf90_targets("${COMMON_RAW_SOURCES}" "${INCLUDES}" "${CPPDEFS}" "" "")
target_sources(common PRIVATE ${SOURCES})

# Build cores!
foreach(CORE IN LISTS CORES)
  build_core(${CORE})
endforeach()
