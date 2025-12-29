#==============================================================================
# Standalone Build Configuration for Emulator Components
#
# Builds all dependencies from E3SM externals:
# - SCORPIO (Parallel I/O)
# - yaml-cpp (from EKAT)
# - Catch2 (from EKAT, for tests)
#
# Requirements: C, C++, Fortran compilers + MPI
#==============================================================================

message(STATUS "")
message(STATUS "=== Emulator Standalone Build Configuration ===")

# Get E3SM root (two levels up from emulator_comps)
get_filename_component(E3SM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(EXTERNALS_DIR "${E3SM_ROOT}/externals")

message(STATUS "E3SM root: ${E3SM_ROOT}")
message(STATUS "Externals: ${EXTERNALS_DIR}")

#------------------------------------------------------------------------------
# Find required tools
#------------------------------------------------------------------------------
find_package(MPI REQUIRED)
message(STATUS "MPI found: ${MPI_CXX_COMPILER}")

#------------------------------------------------------------------------------
# C++17 required
#------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#------------------------------------------------------------------------------
# SCORPIO (Parallel I/O)
# 
# SCORPIO requires NetCDF-C. PnetCDF is optional but recommended.
# On Perlmutter: module load cray-hdf5 cray-netcdf cray-parallel-netcdf
#------------------------------------------------------------------------------
set(SCORPIO_SOURCE_DIR "${EXTERNALS_DIR}/scorpio")

if(EXISTS "${SCORPIO_SOURCE_DIR}/CMakeLists.txt")
  message(STATUS "")
  message(STATUS "--- Building SCORPIO ---")
  message(STATUS "Source: ${SCORPIO_SOURCE_DIR}")
  
  # Look for NetCDF from environment variables (set by Cray modules)
  # Priority: CMake variable > environment variable
  if(NOT NetCDF_C_PATH)
    if(DEFINED ENV{NETCDF_DIR})
      set(NetCDF_C_PATH "$ENV{NETCDF_DIR}" CACHE PATH "Path to NetCDF-C installation")
      message(STATUS "NetCDF from NETCDF_DIR: ${NetCDF_C_PATH}")
    elseif(DEFINED ENV{NETCDF_PATH})
      set(NetCDF_C_PATH "$ENV{NETCDF_PATH}" CACHE PATH "Path to NetCDF-C installation")
      message(STATUS "NetCDF from NETCDF_PATH: ${NetCDF_C_PATH}")
    endif()
  endif()
  
  # Look for PnetCDF from environment variables
  if(NOT PnetCDF_C_PATH)
    if(DEFINED ENV{PNETCDF_DIR})
      set(PnetCDF_C_PATH "$ENV{PNETCDF_DIR}" CACHE PATH "Path to PnetCDF installation")
      message(STATUS "PnetCDF from PNETCDF_DIR: ${PnetCDF_C_PATH}")
    elseif(DEFINED ENV{PNETCDF_PATH})
      set(PnetCDF_C_PATH "$ENV{PNETCDF_PATH}" CACHE PATH "Path to PnetCDF installation")
      message(STATUS "PnetCDF from PNETCDF_PATH: ${PnetCDF_C_PATH}")
    endif()
  endif()
  
  # Check if NetCDF was found
  if(NOT NetCDF_C_PATH)
    message(FATAL_ERROR 
      "NetCDF-C not found! SCORPIO requires NetCDF-C.\n"
      "\n"
      "On Perlmutter (NERSC), load the required modules:\n"
      "  module load cray-hdf5 cray-netcdf cray-parallel-netcdf\n"
      "\n"
      "Then re-run cmake, or set NetCDF_C_PATH explicitly:\n"
      "  cmake .. -DNetCDF_C_PATH=/path/to/netcdf\n")
  endif()
  
  # Set PnetCDF option based on availability
  if(PnetCDF_C_PATH)
    set(WITH_PNETCDF ON CACHE BOOL "Enable PnetCDF" FORCE)
    message(STATUS "PnetCDF enabled: ${PnetCDF_C_PATH}")
  else()
    set(WITH_PNETCDF OFF CACHE BOOL "Disable PnetCDF" FORCE)
    message(STATUS "PnetCDF not found (optional, continuing without it)")
  endif()
  
  # SCORPIO build options - minimal configuration
  set(PIO_ENABLE_FORTRAN ON CACHE BOOL "Enable Fortran in PIO" FORCE)
  set(PIO_ENABLE_TIMING OFF CACHE BOOL "Disable timing" FORCE)
  set(PIO_ENABLE_TESTS OFF CACHE BOOL "Disable PIO tests" FORCE)
  set(PIO_ENABLE_EXAMPLES OFF CACHE BOOL "Disable PIO examples" FORCE)
  set(PIO_ENABLE_TOOLS OFF CACHE BOOL "Disable PIO tools" FORCE)
  set(PIO_USE_MALLOC ON CACHE BOOL "Use malloc" FORCE)
  
  # Add SCORPIO cmake module path so it can find NetCDF/PnetCDF
  list(APPEND CMAKE_MODULE_PATH "${SCORPIO_SOURCE_DIR}/cmake")
  
  # Suppress warnings from SCORPIO (external code we don't control)
  # - Wno-write-strings: bget.cpp has string literal to char* conversions
  set(SCORPIO_COMPILE_FLAGS "-Wno-write-strings")
  
  # Add SCORPIO as subdirectory
  add_subdirectory("${SCORPIO_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/externals/scorpio" EXCLUDE_FROM_ALL)
  
  # Apply warning suppression to SCORPIO target
  if(TARGET pioc)
    target_compile_options(pioc PRIVATE ${SCORPIO_COMPILE_FLAGS})
  endif()
  
  # Set include paths for emulator code
  set(EMULATOR_SCORPIO_INCLUDE_DIRS
    "${SCORPIO_SOURCE_DIR}/src/clib"
    "${CMAKE_BINARY_DIR}/externals/scorpio/src/clib"
    CACHE INTERNAL "SCORPIO include directories")
  
  # Mark SCORPIO as available
  set(EMULATOR_HAS_SCORPIO ON CACHE BOOL "SCORPIO available" FORCE)
  
  message(STATUS "SCORPIO configured successfully")
else()
  message(FATAL_ERROR "SCORPIO not found at ${SCORPIO_SOURCE_DIR}\n"
                      "Ensure E3SM externals are checked out.")
endif()

#------------------------------------------------------------------------------
# yaml-cpp from EKAT
#------------------------------------------------------------------------------
set(YAML_CPP_SOURCE_DIR "${EXTERNALS_DIR}/ekat/extern/yaml-cpp")

if(EXISTS "${YAML_CPP_SOURCE_DIR}/CMakeLists.txt")
  message(STATUS "")
  message(STATUS "--- Building yaml-cpp ---")
  message(STATUS "Source: ${YAML_CPP_SOURCE_DIR}")
  
  # yaml-cpp build options
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Disable yaml-cpp tests" FORCE)
  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Disable yaml-cpp tools" FORCE)
  set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Disable yaml-cpp contrib" FORCE)
  set(YAML_CPP_INSTALL OFF CACHE BOOL "Disable yaml-cpp install" FORCE)
  
  # Suppress warnings from yaml-cpp (external code we don't control)
  # - Wno-deprecated-declarations: std::iterator is deprecated in C++17
  set(YAML_CPP_COMPILE_FLAGS "-Wno-deprecated-declarations")
  
  add_subdirectory("${YAML_CPP_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/externals/yaml-cpp" EXCLUDE_FROM_ALL)
  
  # Apply warning suppression to yaml-cpp target
  if(TARGET yaml-cpp)
    target_compile_options(yaml-cpp PRIVATE ${YAML_CPP_COMPILE_FLAGS})
  endif()
  
  # Set target for linking
  set(EMULATOR_YAML_CPP_TARGET yaml-cpp CACHE INTERNAL "yaml-cpp target")
  
  message(STATUS "yaml-cpp configured successfully")
else()
  message(FATAL_ERROR "yaml-cpp not found in EKAT at ${YAML_CPP_SOURCE_DIR}\n"
                      "Ensure EKAT submodule is initialized:\n"
                      "  git submodule update --init externals/ekat")
endif()

#------------------------------------------------------------------------------
# Catch2 from EKAT (for tests)
#------------------------------------------------------------------------------
if(EMULATOR_BUILD_TESTS)
  set(CATCH2_SOURCE_DIR "${EXTERNALS_DIR}/ekat/extern/Catch2")
  
  if(EXISTS "${CATCH2_SOURCE_DIR}/single_include/catch2/catch.hpp")
    message(STATUS "")
    message(STATUS "--- Catch2 Test Framework ---")
    message(STATUS "Source: ${CATCH2_SOURCE_DIR}")
    
    set(EMULATOR_CATCH2_INCLUDE_DIR "${CATCH2_SOURCE_DIR}/single_include" 
        CACHE INTERNAL "Catch2 include directory")
    
    message(STATUS "Catch2 configured successfully")
  else()
    message(FATAL_ERROR "Catch2 not found in EKAT at ${CATCH2_SOURCE_DIR}\n"
                        "Ensure EKAT submodule is initialized:\n"
                        "  git submodule update --init --recursive externals/ekat")
  endif()
endif()

message(STATUS "")
message(STATUS "=== Standalone Configuration Complete ===")
message(STATUS "")
