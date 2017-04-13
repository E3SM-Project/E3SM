cmake_minimum_required (VERSION 2.8.12)
include (CheckFunctionExists)
project (PIOC C)

#==============================================================================
#  DEFINE THE TARGET
#==============================================================================

add_library (pioc topology.c pio_file.c pioc_support.c pio_lists.c
  pioc.c pioc_sc.c pio_spmd.c pio_rearrange.c pio_nc4.c bget.c
  pio_nc.c pio_put_nc.c pio_get_nc.c pio_getput_int.c pio_msg.c pio_varm.c
  pio_darray.c pio_darray_int.c)

# set up include-directories
include_directories(
  "${PROJECT_SOURCE_DIR}"   # to find foo/foo.h
  "${PROJECT_BINARY_DIR}")  # to find foo/config.h

# Include the clib source directory
target_include_directories (pioc
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# System and compiler CPP directives
target_compile_definitions (pioc
  PUBLIC ${CMAKE_SYSTEM_DIRECTIVE})
target_compile_definitions (pioc
  PUBLIC ${CMAKE_C_COMPILER_DIRECTIVE})

# Compiler-specific compiler options
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options (pioc
    PRIVATE -std=c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  target_compile_options (pioc
    PRIVATE -c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  target_compile_options (pioc
    PRIVATE -std=c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options (pioc
    PRIVATE -std=c99)
endif()

#==============================================================================
#  DEFINE THE INSTALL
#==============================================================================

# Library
install (TARGETS pioc DESTINATION lib)

# Include/Header File
install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/pio.h DESTINATION include)

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================

#===== MPI =====
if (PIO_USE_MPISERIAL)
  find_package (MPISERIAL COMPONENTS C REQUIRED)
  if (MPISERIAL_C_FOUND)
    target_compile_definitions (pioc
      PRIVATE MPI_SERIAL)
    target_include_directories (pioc
      PUBLIC ${MPISERIAL_C_INCLUDE_DIRS})
    target_link_libraries (pioc
      PUBLIC ${MPISERIAL_C_LIBRARIES})

    set (WITH_PNETCDF FALSE)
  endif ()
else ()
  find_package (MPI REQUIRED)
endif ()

#===== GPTL =====
if (PIO_ENABLE_TIMING)
  find_package (GPTL COMPONENTS C QUIET)
  if (GPTL_C_FOUND)
    message (STATUS "Found GPTL C: ${GPTL_C_LIBRARIES}")
    target_include_directories (pioc
      PUBLIC ${GPTL_C_INCLUDE_DIRS})
    target_link_libraries (pioc
      PUBLIC ${GPTL_C_LIBRARIES})
  else ()
    message (STATUS "Using internal GPTL C library for timing")
    target_link_libraries (pioc
      PUBLIC gptl)
  endif ()
  target_compile_definitions (pioc
    PUBLIC TIMING)
endif ()

#===== NetCDF-C =====
find_package (NetCDF "4.3.3" COMPONENTS C)
if (NetCDF_C_FOUND)
  target_include_directories (pioc
    PUBLIC ${NetCDF_C_INCLUDE_DIRS})
  target_compile_definitions (pioc
    PUBLIC _NETCDF)
  target_link_libraries (pioc
    PUBLIC ${NetCDF_C_LIBRARIES})
  if (${NetCDF_C_HAS_PARALLEL})
    target_compile_definitions (pioc
      PUBLIC _NETCDF4)
  endif ()
  if (${NetCDF_C_LOGGING_ENABLED})
    target_compile_definitions (pioc
      PUBLIC NETCDF_C_LOGGING_ENABLED)
    # netcdf.h needs this to be defined to use netCDF logging.
    target_compile_definitions (pioc
      PUBLIC LOGGING)
  endif()
else ()
  target_compile_definitions (pioc
    PUBLIC _NONETCDF)
endif ()

#===== PnetCDF-C =====
if (WITH_PNETCDF)
  find_package (PnetCDF "1.7.0" COMPONENTS C)
endif ()
if (PnetCDF_C_FOUND)
  target_include_directories (pioc
    PUBLIC ${PnetCDF_C_INCLUDE_DIRS})
  target_compile_definitions (pioc
    PUBLIC _PNETCDF)
  target_link_libraries (pioc
    PUBLIC ${PnetCDF_C_LIBRARIES})

  # Check library for varn functions
  set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_C_LIBRARY})
  check_function_exists (ncmpi_get_varn PnetCDF_C_HAS_VARN)
  if (PnetCDF_C_HAS_VARN)
    target_compile_definitions(pioc
      PUBLIC USE_PNETCDF_VARN
      PUBLIC USE_PNETCDF_VARN_ON_READ)
  endif()
else ()
  target_compile_definitions (pioc
    PUBLIC _NOPNETCDF)
endif ()

#===== Add EXTRAs =====
target_include_directories (pioc
  PUBLIC ${PIO_C_EXTRA_INCLUDE_DIRS})
target_link_libraries (pioc
  PUBLIC ${PIO_C_EXTRA_LIBRARIES})
target_compile_options (pioc
  PRIVATE ${PIO_C_EXTRA_COMPILE_OPTIONS})
target_compile_definitions (pioc
  PUBLIC ${PIO_C_EXTRA_COMPILE_DEFINITIONS})
if (PIO_C_EXTRA_LINK_FLAGS)
  set_target_properties(pioc PROPERTIES
    LINK_FLAGS ${PIO_C_EXTRA_LINK_FLAGS})
endif ()

#===== Check for necessities =====
if (NOT PnetCDF_C_FOUND AND NOT NetCDF_C_FOUND)
  message (FATAL_ERROR "Must have PnetCDF and/or NetCDF C libraries")
endif ()

include(CheckTypeSize)
check_type_size("size_t" SIZEOF_SIZE_T)
CHECK_TYPE_SIZE("long long" SIZEOF_LONG_LONG)
if (NOT ${SIZEOF_SIZE_T} EQUAL ${SIZEOF_LONG_LONG})
  message (FATAL_ERROR "size_t and long long must be the same size!")
endif ()
