#
# This is the shared entry point for building MPAS components with CMake.
#
# MPAS_BUILD_MODE selects the mode-specific pieces:
#  * E3SM : coupled builds through CIME (set by build_mpas_models() in
#           components/cmake/build_mpas_model.cmake). This is the default.
#  * STANDALONE : standalone builds (not yet implemented)
#
# The following CMake variables are expected to be defined:
#  * CORES : A list of CORES to build, comma-separated (e.g. "ocean,seaice,landice")
#  * In E3SM mode, whatever CIME settings are setting to correctly resolve the
#    ${CASEROOT}/Macros.cmake file
#    - COMPILER, DEBUG, MPILIB, MACH, OS
#

if (NOT MPAS_BUILD_MODE)
  set(MPAS_BUILD_MODE "E3SM")
endif()
if (MPAS_BUILD_MODE STREQUAL "STANDALONE")
  message(FATAL_ERROR "MPAS_BUILD_MODE=STANDALONE is not yet implemented")
elseif (NOT MPAS_BUILD_MODE STREQUAL "E3SM")
  message(FATAL_ERROR "Unrecognized MPAS_BUILD_MODE: ${MPAS_BUILD_MODE} (expected E3SM or STANDALONE)")
endif()

if (MPAS_BUILD_MODE STREQUAL "E3SM")
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
endif()

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
if (MPAS_BUILD_MODE STREQUAL "E3SM")
  list(APPEND CPPDEFS "-DMPAS_NO_LOG_REDIRECT" "-DMPAS_NO_ESMF_INIT" "-DMPAS_ESM_SHR_CONST" "-DMPAS_PERF_MOD_TIMERS")
endif()
list(APPEND CPPDEFS "${MODEL_FORMULATION}" "${FILE_OFFSET}" "${ZOLTAN_DEFINE}" "-D_MPI" "-DMPAS_NAMELIST_SUFFIX=${NAMELIST_SUFFIX}" "-DMPAS_EXE_NAME=${EXE_NAME}")
if (PIO_VERSION STREQUAL "2")
  list(APPEND CPPDEFS "-DUSE_PIO2")
  if (TARGET adios2::adios2)
    list(APPEND CPPDEFS "-DHAVE_ADIOS")
  endif()
  if (TARGET hdf5)
    list(APPEND CPPDEFS "-DHAVE_HDF5")
  endif()
endif()

if (DEBUG)
  list(APPEND CPPDEFS "-DMPAS_DEBUG")
endif()
if (compile_threaded)
  list(APPEND CPPDEFS "-DMPAS_OPENMP")
endif()

set(INCLUDES)
if (MPAS_BUILD_MODE STREQUAL "E3SM")
  list(APPEND INCLUDES "${INSTALL_SHAREDPATH}/include")
endif()
list(APPEND INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/ezxml" "${CMAKE_BINARY_DIR}/framework" "${CMAKE_BINARY_DIR}/operators")

if (MPAS_BUILD_MODE STREQUAL "E3SM" AND USE_KOKKOS)
  include(${INSTALL_SHAREDPATH}/kokkos_generated_settings.cmake)
  string (REPLACE ";" " " KOKKOS_CXXFLAGS_STR "${KOKKOS_CXXFLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${KOKKOS_CXXFLAGS_STR}")
endif()

# Track the files included by preprocessed sources and registries. DEPFILE
# needs CMake 3.20 or later with Makefile generators. Older versions fall back
# to depending only on the source file itself.
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.20)
  set(MPAS_USE_DEPFILE TRUE)
  # Avoids a policy warning with Ninja. Depfile paths are absolute, so the
  # OLD and NEW behaviors give the same result.
  cmake_policy(SET CMP0116 NEW)
else()
  set(MPAS_USE_DEPFILE FALSE)
endif()

# Include custom cmake libraries used for mpas
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_utils.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_core.cmake)

# Add tools
add_subdirectory(tools)

# Gather sources that are needed for all cores into "common" library

set(COMMON_RAW_SOURCES external/ezxml/ezxml.c)
include(${CMAKE_CURRENT_SOURCE_DIR}/framework/framework.cmake)
if (MPAS_BUILD_MODE STREQUAL "E3SM" AND COMP_INTERFACE STREQUAL "moab")
  list(APPEND COMMON_RAW_SOURCES framework/mpas_moabmesh.F)
  list(APPEND INCLUDES ${MOAB_INCLUDE_DIRS})
  list(APPEND CPPDEFS "-DHAVE_MOAB")
endif()
include(${CMAKE_CURRENT_SOURCE_DIR}/operators/operators.cmake)

add_library(common)
target_compile_definitions(common PRIVATE ${CPPDEFS})
target_include_directories(common PRIVATE ${INCLUDES})
if (MPAS_BUILD_MODE STREQUAL "E3SM")
  target_link_libraries(common PUBLIC spio csm_share)
endif()

genf90_targets("${COMMON_RAW_SOURCES}" "${INCLUDES}" "${CPPDEFS}" "" "")
target_sources(common PRIVATE ${SOURCES})

# Build cores!
foreach(CORE IN LISTS CORES)
  build_core(${CORE})
endforeach()
