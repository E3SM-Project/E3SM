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

include_directories (${CMAKE_CURRENT_SOURCE_DIR})

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set_target_properties (pioc PROPERTIES COMPILE_FLAGS -std=c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    set_target_properties (pioc PROPERTIES COMPILE_FLAGS -c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set_target_properties (pioc PROPERTIES COMPILE_FLAGS -std=c99)
endif()

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================
    
#===== GPTL =====
if (PIO_ENABLE_TIMING)
    add_definitions(-DTIMING)
    include_directories(${GPTL_INCLUDE_DIRS})
    target_link_libraries (pioc gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)
    
#===== NetCDF-C =====
find_package (NetCDF "4.3.3" COMPONENTS C)
if (NetCDF_C_FOUND)
    include_directories (${NetCDF_C_INCLUDE_DIRS})
    add_definitions (-D_NETCDF)
    target_link_libraries (pioc ${NetCDF_C_LIBRARIES})
    if (${NetCDF_C_HAS_PARALLEL})
        add_definitions (-D_NETCDF4)
    endif ()
else ()
    add_definitions (-D_NONETCDF)
endif ()

#===== PnetCDF-C =====
find_package (PnetCDF "1.6" COMPONENTS C)
if (PnetCDF_C_FOUND)
    include_directories (${PnetCDF_C_INCLUDE_DIRS})
    add_definitions (-D_PNETCDF)
    target_link_libraries (pioc ${PnetCDF_C_LIBRARIES})

    # Check library for varn functions
    include(CheckFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_C_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
        add_definitions(-DUSE_PNETCDF_VARN -DUSE_PNETCDF_VARN_ON_READ)
    endif()  
else ()
    add_definitions (-D_NOPNETCDF)
endif ()

#===== Add EXTRAs =====
include_directories (${PIO_C_EXTRA_INCLUDE_DIRS})
target_link_libraries (pioc ${PIO_C_EXTRA_LIBRARIES})
add_definitions (${PIO_C_EXTRA_DEFINITIONS})
if (PIO_C_EXTRA_LINK_FLAGS)
    set_target_properties(pioc PROPERTIES 
        LINK_FLAGS ${PIO_C_EXTRA_LINK_FLAGS})
endif ()

#===== Check for necessities =====
if (NOT PnetCDF_C_FOUND AND NOT NetCDF_C_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF C libraries")
endif ()

#==============================================================================
#  DEFINE PARAMETERS FOR TARGETS DEPENDENT ON THIS
#==============================================================================

get_target_property(PIO_C_INCLUDE_DIRS pioc INCLUDE_DIRECTORIES)
message ("PIO_C_INCLUDE_DIRS = ${PIO_C_INCLUDE_DIRS}")
set (PIO_C_INCLUDE_DIRS ${PIO_C_INCLUDE_DIRS}
     CACHE STRING "PIO C Include Directory")
