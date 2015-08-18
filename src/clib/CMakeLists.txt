cmake_minimum_required (VERSION 2.8.10)
project (PIOC C)

#==============================================================================
#  DEFINE THE TARGET
#==============================================================================

set (PIO_C_SRCS topology.c
                pio_nc.c
                pio_file.c
                pioc_support.c
                pio_lists.c
                pioc.c
                pioc_sc.c
                pio_spmd.c
                pio_rearrange.c
                pio_darray.c
                pio_put_nc.c 
                pio_get_nc.c
                bget.c)

add_library (pioc ${PIO_C_SRCS})

set (PIO_C_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
set (PIO_C_DEFINITIONS)
set (PIO_C_COMPILE_FLAGS)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    list (APPEND PIO_C_COMPILE_FLAGS -std=c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    list (APPEND PIO_C_COMPILE_FLAGS -std=c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    list (APPEND PIO_C_COMPILE_FLAGS -c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    list (APPEND PIO_C_COMPILE_FLAGS -std=c99)
endif()

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================
    
#===== GPTL =====
if (PIO_ENABLE_TIMING)
    list (APPEND PIO_C_INCLUDE_DIRS ${GPTL_INCLUDE_DIRS})
    list (APPEND PIO_C_DEFINITIONS TIMING ${GPTL_DEFINITIONS})
    target_link_libraries (pioc gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)
    
#===== NetCDF-C =====
find_package (NetCDF "4.3.3" COMPONENTS C)
if (NetCDF_C_FOUND)
    list (APPEND PIO_C_INCLUDE_DIRS ${NetCDF_C_INCLUDE_DIRS})
    list (APPEND PIO_C_DEFINITIONS _NETCDF)
    if (NetCDF_C_HAS_PARALLEL)
        list (APPEND PIO_C_DEFINITIONS _NETCDF4)
    endif ()
    target_link_libraries (pioc ${NetCDF_C_LIBRARIES})
else ()
    list (APPEND PIO_C_DEFINITIONS _NONETCDF)
endif ()

#===== PnetCDF-C =====
find_package (PnetCDF "1.6" COMPONENTS C)
if (PnetCDF_C_FOUND)
    list (APPEND PIO_C_INCLUDE_DIRS ${PnetCDF_C_INCLUDE_DIRS})
    list (APPEND PIO_C_DEFINITIONS _PNETCDF)
    target_link_libraries (pioc ${PnetCDF_C_LIBRARIES})

    # Check library for varn functions
    include(CheckFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_C_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
        list (APPEND PIO_C_DEFINITIONS USE_PNETCDF_VARN USE_PNETCDF_VARN_ON_READ)
    endif()  
else ()
    list (APPEND PIO_C_DEFINITIONS _NOPNETCDF)
endif ()

#===== Add EXTRAs =====
list (APPEND PIO_C_INCLUDE_DIRS ${PIO_C_EXTRA_INCLUDE_DIRS})
list (APPEND PIO_C_DEFINITIONS ${PIO_C_EXTRA_DEFINITIONS})
list (APPEND PIO_C_COMPILE_FLAGS ${PIO_C_EXTRA_COMPILE_FLAGS})
target_link_libraries (pioc ${PIO_C_EXTRA_LIBRARIES})

#===== Check for necessities =====
if (NOT PnetCDF_C_FOUND AND NOT NetCDF_C_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF C libraries")
endif ()

#==============================================================================
#  SET TARGET PROPERTIES
#==============================================================================

set_target_properties (pioc PROPERTIES
    INCLUDE_DIRECTORIES "${PIO_C_INCLUDE_DIRS}"
    COMPILE_DEFINITIONS "${PIO_C_DEFINITIONS}"
    COMPILE_FLAGS "${PIO_C_COMPILE_FLAGS}")

# Cache for out-of-scope use
set (PIO_C_INCLUDE_DIRS "${PIO_C_INCLUDE_DIRS}"
     CACHE STRING "PIO C include directories")
set (PIO_C_DEFINITIONS "${PIO_C_DEFINITIONS}"
     CACHE STRING "PIO C compiler definitions")
set (PIO_C_COMPILE_FLAGS "${PIO_C_COMPILE_FLAGS}"
     CACHE STRING "PIO C compiler flags")

