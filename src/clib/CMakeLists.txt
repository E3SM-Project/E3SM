cmake_minimum_required (VERSION 3.0.0)
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

target_include_directories (pioc
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options (pioc
        PRIVATE -std=c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    target_compile_options (pioc
        PRIVATE -c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    target_compile_options (pioc
        PRIVATE -std=c99)
endif()

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================
    
#===== GPTL =====
if (ENABLE_TIMING)
    target_compile_definitions (pioc
        PRIVATE TIMING)
    target_link_libraries (pioc
        PRIVATE gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)
#target_include_directories (pioc 
#    PUBLIC ${MPI_C_INCLUDE_PATH})
target_compile_definitions (pioc 
    PRIVATE USEMPIIO)
#target_link_libraries (pioc
#    PRIVATE ${MPI_C_LIBRARIES})

#===== NetCDF-C =====
find_package (NetCDF_C)
if (NetCDF_C_FOUND)
    target_include_directories (pioc 
        PUBLIC ${NetCDF_C_INCLUDE_DIRS})
    target_compile_definitions (pioc 
        PRIVATE _NETCDF)
    target_link_libraries (pioc
        PRIVATE ${NetCDF_C_LIBRARIES})
    if (EXISTS ${NetCDF_C_INCLUDE_DIR}/netcdf_par.h)
        target_compile_definitions (pioc 
            PRIVATE _NETCDF4)
    endif ()
else ()
    target_compile_definitions (pioc 
        PRIVATE _NONETCDF)
endif ()

#===== PnetCDF =====
find_package (PnetCDF_C)
if (PnetCDF_C_FOUND)
    target_include_directories (pioc 
        PUBLIC ${PnetCDF_C_INCLUDE_DIRS})
    target_compile_definitions (pioc 
        PRIVATE _PNETCDF)
    target_link_libraries (pioc
        PRIVATE ${PnetCDF_C_LIBRARIES})
        
    # Check library for varn functions
    include(CheckFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_C_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
        target_compile_definitions(pioc
            PRIVATE USE_PNETCDF_VARN
            PRIVATE USE_PNETCDF_VARN_ON_READ)
    endif()  
else ()
    target_compile_definitions (pioc 
        PRIVATE _NOPNETCDF)
endif ()

#===== Check for necessities =====
if (NOT PnetCDF_C_FOUND AND NOT NetCDF_C_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF C libraries")
endif ()
