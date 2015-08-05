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
    target_link_libraries (pioc gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)
target_include_directories (pioc 
    PUBLIC ${MPI_C_INCLUDE_PATH})
target_compile_definitions (pioc 
    PRIVATE USEMPIIO)
target_link_libraries (pioc
    ${MPI_C_LIBRARIES})

#===== NetCDF-C =====
find_package (NetCDF)
if (NetCDF_C_FOUND)
    target_include_directories (pioc 
        PUBLIC ${NetCDF_C_INCLUDE_DIRS})
    target_compile_definitions (pioc 
        PRIVATE _NETCDF)
    target_link_libraries (pioc
        ${NetCDF_C_LIBRARIES})
    if (EXISTS ${NetCDF_C_INCLUDE_DIR}/netcdf_par.h)
        target_compile_definitions (pioc 
            PRIVATE _NETCDF4)
    endif ()
else ()
    target_compile_definitions (pioc 
        PRIVATE _NONETCDF)
endif ()

#===== PnetCDF =====
find_package (PnetCDF)
target_include_directories (pioc 
    PUBLIC ${PnetCDF_C_INCLUDE_DIRS})
target_compile_definitions (pioc 
    PRIVATE _PNETCDF
    PRIVATE USE_PNETCDF_VARN
    PRIVATE USE_PNETCDF_VARN_ON_READ)
target_link_libraries (pioc
    ${PnetCDF_C_LIBRARIES})
