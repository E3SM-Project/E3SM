cmake_minimum_required (VERSION 2.8.12)
project (PIOC C)

#==============================================================================
#  DEFINE THE TARGET
#==============================================================================

set (PIO_C_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/topology.c
                ${CMAKE_CURRENT_BINARY_DIR}/pio_nc.c
                ${CMAKE_CURRENT_SOURCE_DIR}/pio_file.c
                ${CMAKE_CURRENT_SOURCE_DIR}/pioc_support.c
                ${CMAKE_CURRENT_SOURCE_DIR}/pio_lists.c
                ${CMAKE_CURRENT_SOURCE_DIR}/pioc.c
                ${CMAKE_CURRENT_SOURCE_DIR}/pioc_sc.c
                ${CMAKE_CURRENT_SOURCE_DIR}/pio_spmd.c
                ${CMAKE_CURRENT_SOURCE_DIR}/pio_rearrange.c
                ${CMAKE_CURRENT_SOURCE_DIR}/pio_darray.c
                ${CMAKE_CURRENT_BINARY_DIR}/pio_put_nc.c 
                ${CMAKE_CURRENT_BINARY_DIR}/pio_get_nc.c
                ${CMAKE_CURRENT_SOURCE_DIR}/bget.c)

add_library (pioc ${PIO_C_SRCS})

target_include_directories (pioc
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

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
#  DEFINE THE DEPENDENCIES
#==============================================================================
    
#===== GPTL =====
if (PIO_ENABLE_TIMING)
    target_compile_definitions (pioc
        PUBLIC TIMING)
    target_link_libraries (pioc
        PUBLIC gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)
    
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
else ()
    target_compile_definitions (pioc 
        PUBLIC _NONETCDF)
endif ()

#===== PnetCDF-C =====
find_package (PnetCDF "1.6" COMPONENTS C)
if (PnetCDF_C_FOUND)
    target_include_directories (pioc 
        PUBLIC ${PnetCDF_C_INCLUDE_DIRS})
    target_compile_definitions (pioc 
        PUBLIC _PNETCDF)
    target_link_libraries (pioc
        PUBLIC ${PnetCDF_C_LIBRARIES})

    # Check library for varn functions
    include(CheckFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_C_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
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

#===== Generate missing sources from templates =====
# If both NetCDF and PnetCDF were found, then create templates
if (PnetCDF_C_FOUND AND NetCDF_C_FOUND AND PIO_GENERATE_TEMPLATE_SOURCES)

    # Find perl
    find_package (Perl REQUIRED)
    if (Perl_FOUND)
        execute_process (
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pio_c_put_template.c
                                             ${CMAKE_CURRENT_BINARY_DIR}/pio_c_put_template.c
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        execute_process (
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pio_c_get_template.c
                                             ${CMAKE_CURRENT_BINARY_DIR}/pio_c_get_template.c
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        execute_process (
            COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ncputgetparser.pl 
                                       ${NetCDF_C_INCLUDE_DIR} ${PnetCDF_C_INCLUDE_DIR}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        execute_process (
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pio_c_template.c
                                             ${CMAKE_CURRENT_BINARY_DIR}/pio_c_template.c
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        execute_process (
            COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ncparser.pl
                                       ${NetCDF_C_INCLUDE_DIR} ${PnetCDF_C_INCLUDE_DIR}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    else ()
        message (FATAL_ERROR "Need Perl to create PIO C source files from templates")
    endif ()
    
# If both NetCDF and PnetCDF are NOT found, then just copy existing source files
else ()
        execute_process (
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pio_put_nc.c
                                             ${CMAKE_CURRENT_BINARY_DIR}/pio_put_nc.c
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pio_get_nc.c
                                             ${CMAKE_CURRENT_BINARY_DIR}/pio_get_nc.c
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pio_nc.c
                                             ${CMAKE_CURRENT_BINARY_DIR}/pio_nc.c
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif ()


#===== Check for necessities =====
if (NOT PnetCDF_C_FOUND AND NOT NetCDF_C_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF C libraries")
endif ()
