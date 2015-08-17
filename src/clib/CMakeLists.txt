cmake_minimum_required (VERSION 2.8)
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
set_target_properties(pioc PROPERTIES 
                      INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})
#target_include_directories (pioc
#    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set_target_properties(pioc PROPERTIES 
                          COMPILE_FLAGS -std=c99)
#    target_compile_options (pioc
#        PRIVATE -std=c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    set_target_properties(pioc PROPERTIES 
                          COMPILE_FLAGS -c99)
#    target_compile_options (pioc
#        PRIVATE -c99)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set_target_properties(pioc PROPERTIES 
                          COMPILE_FLAGS -std=c99)
#    target_compile_options (pioc
#        PRIVATE -std=c99)
endif()

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================
    
#===== GPTL =====
if (PIO_ENABLE_TIMING)
#    target_compile_definitions (pioc
#        PRIVATE TIMING)
    set_target_properties(pioc PROPERTIES 
                          COMPILE_DEFINITIONS -std=c99)
    target_link_libraries (pioc
        PUBLIC gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)
    
#===== NetCDF-C =====
find_package (NetCDF "4.3.3" COMPONENTS C)
if (NetCDF_C_FOUND)
    set_target_properties(pioc PROPERTIES 
                          INCLUDE_DIRECTORIES ${NetCDF_C_INCLUDE_DIRS})
#    target_include_directories (pioc 
#        PUBLIC ${NetCDF_C_INCLUDE_DIRS})
    set_target_properties(pioc PROPERTIES 
                          COMPILE_DEFINITIONS _NETCDF)
#    target_compile_definitions (pioc 
#        PRIVATE _NETCDF)
    target_link_libraries (pioc
        PUBLIC ${NetCDF_C_LIBRARIES})
    if (${NetCDF_C_HAS_PARALLEL})
        set_target_properties(pioc PROPERTIES 
                              COMPILE_DEFINITIONS _NETCDF4)
#        target_compile_definitions (pioc 
#            PRIVATE _NETCDF4)
    endif ()
else ()
    set_target_properties(pioc PROPERTIES 
                          COMPILE_DEFINITIONS _NONETCDF)
#    target_compile_definitions (pioc 
#        PRIVATE _NONETCDF)
endif ()

#===== PnetCDF-C =====
find_package (PnetCDF "1.6" COMPONENTS C)
if (PnetCDF_C_FOUND)
    set_target_properties(pioc PROPERTIES 
                          INCLUDE_DIRECTORIES ${PnetCDF_C_INCLUDE_DIRS})
#    target_include_directories (pioc 
#        PUBLIC ${PnetCDF_C_INCLUDE_DIRS})
    set_target_properties(pioc PROPERTIES 
                          COMPILE_DEFINITIONS _PNETCDF)
#    target_compile_definitions (pioc 
#        PRIVATE _PNETCDF)
    target_link_libraries (pioc
        PUBLIC ${PnetCDF_C_LIBRARIES})

    # Check library for varn functions
    include(CheckFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_C_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
        set_target_properties(pioc PROPERTIES 
                              COMPILE_DEFINITIONS USE_PNETCDF_VARN USE_PNETCDF_VARN_ON_READ)
#        target_compile_definitions(pioc
#            PRIVATE USE_PNETCDF_VARN
#            PRIVATE USE_PNETCDF_VARN_ON_READ)
    endif()  
else ()
    set_target_properties(pioc PROPERTIES 
                          COMPILE_DEFINITIONS _NOPNETCDF)
#    target_compile_definitions (pioc 
#        PRIVATE _NOPNETCDF)
endif ()

#===== Add EXTRAs =====
#target_include_directories (pioc
#    PUBLIC ${PIO_C_EXTRA_INCLUDE_DIRS})
#target_link_libraries (pioc
#    PUBLIC ${PIO_C_EXTRA_LIBRARIES})
#target_compile_options (pioc
#    PUBLIC ${PIO_C_EXTRA_COMPILE_OPTIONS})
#target_compile_definitions (pioc
#    PUBLIC ${PIO_C_EXTRA_COMPILE_DEFINITIONS})
#if (PIO_C_EXTRA_LINK_FLAGS)
#    set_target_properties(pioc PROPERTIES 
#        LINK_FLAGS ${PIO_C_EXTRA_LINK_FLAGS})
#endif ()

#===== Check for necessities =====
if (NOT PnetCDF_C_FOUND AND NOT NetCDF_C_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF C libraries")
endif ()
