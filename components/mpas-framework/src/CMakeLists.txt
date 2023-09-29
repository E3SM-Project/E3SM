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

set(CMAKE_C_COMPILER ${MPICC})
set(CMAKE_CXX_COMPILER ${MPICXX})
set(CMAKE_Fortran_COMPILER ${MPIFC})
set(CMAKE_EXE_LINKER_FLAGS "${LDFLAGS}")
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Set up CPPDEFS
set(FILE_OFFSET "-DOFFSET64BIT")
if (CPPDEFS)
  separate_arguments(CPPDEFS UNIX_COMMAND "${CPPDEFS}")
endif()
list(APPEND CPPDEFS "-DMPAS_NO_LOG_REDIRECT" "-DUSE_PIO2" "-DHAVE_ADIOS" "-DMPAS_NO_ESMF_INIT" "-DMPAS_ESM_SHR_CONST" "-DMPAS_PERF_MOD_TIMERS" "${MODEL_FORMULATION}" "${FILE_OFFSET}" "${ZOLTAN_DEFINE}" "-D_MPI" "-DMPAS_NAMELIST_SUFFIX=${NAMELIST_SUFFIX}" "-DMPAS_EXE_NAME=${EXE_NAME}")
if (DEBUG)
  list(APPEND CPPDEFS "-DMPAS_DEBUG")
endif()
if (compile_threaded)
  list(APPEND CPPDEFS "-DMPAS_OPENMP")
endif()

set(INCLUDES "${INSTALL_SHAREDPATH}/include" "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share" "${INSTALL_SHAREDPATH}/pio" "${CMAKE_CURRENT_SOURCE_DIR}/external/ezxml" "${CMAKE_BINARY_DIR}/framework" "${CMAKE_BINARY_DIR}/operators")

if (USE_KOKKOS)
  include(${INSTALL_SHAREDPATH}/kokkos_generated_settings.cmake)
  string (REPLACE ";" " " KOKKOS_CXXFLAGS_STR "${KOKKOS_CXXFLAGS}")
  set(CXXFLAGS "${CXXFLAGS} ${KOKKOS_CXXFLAGS_STR}")
endif()

set(CMAKE_Fortran_FLAGS "${FFLAGS}")
set(CMAKE_C_FLAGS "${CFLAGS}")
set(CMAKE_CXX_FLAGS "${CXXFLAGS}")

# Include custom cmake libraries used for mpas
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_utils.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_core.cmake)

# Add tools
add_subdirectory(tools)

# Gather sources that are needed for all cores into "common" library

set(COMMON_RAW_SOURCES external/ezxml/ezxml.c)
include(${CMAKE_CURRENT_SOURCE_DIR}/framework/framework.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/operators/operators.cmake)

add_library(common)
target_compile_definitions(common PRIVATE ${CPPDEFS})
target_include_directories(common PRIVATE ${INCLUDES})
find_package(PIO REQUIRED)
target_link_libraries(common PUBLIC spio)

genf90_targets("${COMMON_RAW_SOURCES}" "${INCLUDES}" "${CPPDEFS}" "" "")
target_sources(common PRIVATE ${SOURCES})

# Build cores!
foreach(CORE IN LISTS CORES)
  build_core(${CORE})
endforeach()
