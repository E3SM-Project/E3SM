cmake_minimum_required (VERSION 3.0.0)
project (PIOF Fortran)

#==============================================================================
#  APPEND TO CMAKE MODULE PATH
#==============================================================================

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#==============================================================================
#  DEFINE THE TARGET
#==============================================================================

set (PIO_Fortran_SRCS pio_nf.F90
                      pio.F90
                      pio_kinds.F90
                      pio_types.F90
                      piolib_mod.F90 
                      pio_support.F90)
                      
set (PIO_GenF90_SRCS pionfatt_mod.F90
                     pionfput_mod.F90
                     pionfget_mod.F90
                     piodarray.F90)

foreach (SRC_FILE IN LISTS PIO_GenF90_SRCS)
add_custom_command (OUTPUT ${SRC_FILE}
                    COMMAND ${GENF90_DIR}/genf90.pl
                            ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in > ${SRC_FILE}
                    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in)
endforeach ()

add_library (piof ${PIO_Fortran_SRCS} ${PIO_GenF90_SRCS})
if (NOT ENABLE_FORTRAN)
    set_target_properties(piof PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif ()

target_include_directories (piof
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options (piof
        PRIVATE -ffree-line-length-none)
endif()

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================

#===== GPTL =====
if (ENABLE_TIMING)
    target_compile_definitions (piof
        PRIVATE TIMING)
    target_link_libraries (piof 
        PUBLIC gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)
target_compile_definitions (piof 
    PRIVATE USEMPIIO)

#===== NetCDF-Fortran =====
if ("${MPI_Fortran_LIBRARIES}" MATCHES "libnetcdf")
    message (STATUS "NetCDF_Fortran found in MPI_Fortran wrapper")
    set (NetCDF_Fortran_FOUND TRUE)
    target_compile_definitions (piof 
        PRIVATE _NETCDF)
    find_path (NetCDF_Fortran_IN_MPI_INCLUDES
               NAMES netcdf_par.h
               PATHS ${MPI_Fortran_INCLUDE_PATH})
    if (NetCDF_Fortran_IN_MPI_INCLUDES)
        target_compile_definitions (piof
            PRIVATE _NETCDF4)
    endif ()
else ()
    find_package (NetCDF COMPONENTS Fortran)
    if (NetCDF_Fortran_FOUND)
        target_include_directories (piof 
            PUBLIC ${NetCDF_Fortran_INCLUDE_DIRS})
        target_compile_definitions (piof 
            PRIVATE _NETCDF)
        target_link_libraries (piof
            PUBLIC ${NetCDF_Fortran_LIBRARIES})
        if (EXISTS ${NetCDF_Fortran_INCLUDE_DIR}/netcdf_par.h)
            target_compile_definitions (piof
                PRIVATE _NETCDF4)
        endif ()
    else ()
        target_compile_definitions (piof 
            PRIVATE _NONETCDF)
    endif ()
endif ()

#===== PnetCDF =====
if ("${MPI_Fortran_LIBRARIES}" MATCHES "libpnetcdf")
    message (STATUS "PnetCDF_Fortran found in MPI_Fortran wrapper")
    set (PnetCDF_Fortran_FOUND TRUE)
    target_compile_definitions (piof 
        PRIVATE _PNETCDF)
    set (PnetCDF_Fortran_LIBRARY ${MPI_Fortran_LIBRARIES})
else ()
    find_package (PnetCDF COMPONENTS Fortran)
    if (PnetCDF_Fortran_FOUND)
        target_include_directories (piof 
            PUBLIC ${PnetCDF_Fortran_INCLUDE_DIRS})
        target_compile_definitions (piof 
            PRIVATE _PNETCDF)
        target_link_libraries (piof
            PUBLIC ${PnetCDF_Fortran_LIBRARIES})
    endif ()
endif ()
if (PnetCDF_Fortran_FOUND)
    # Check library for varn functions
    include(CheckFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_Fortran_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
        target_compile_definitions(piof
            PRIVATE USE_PNETCDF_VARN
            PRIVATE USE_PNETCDF_VARN_ON_READ)
    endif()  
else ()
    target_compile_definitions (piof 
        PRIVATE _NOPNETCDF)
endif ()

#===== Check for necessities =====
if (NOT PnetCDF_Fortran_FOUND AND NOT NetCDF_Fortran_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF Fortran libraries")
endif ()
