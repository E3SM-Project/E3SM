cmake_minimum_required (VERSION 3.0.0)
project (PIOF Fortran)

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
        PUBLIC -ffree-line-length-none)
endif()

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================

#===== GPTL =====
if (ENABLE_TIMING)
    target_compile_definitions (piof
        PUBLIC TIMING)
    target_link_libraries (piof gptl)
endif ()

#===== MPI =====
find_package (MPI)
target_include_directories (piof 
    PUBLIC ${MPI_Fortran_INCLUDE_PATH})
target_compile_definitions (piof 
    PUBLIC USEMPIIO)
target_link_libraries (piof
    ${MPI_Fortran_LIBRARIES})

#===== NetCDF =====
target_include_directories (piof 
    PUBLIC $ENV{NETCDF}/include)
target_compile_definitions (piof 
    PUBLIC _NETCDF
    PUBLIC _NETCDF4)
target_link_libraries (piof
    $ENV{NETCDF}/lib/libnetcdff.a)

#===== PnetCDF =====
target_include_directories (piof 
    PUBLIC $ENV{PNETCDF}/include)
target_compile_definitions (piof 
    PUBLIC _PNETCDF
    PUBLIC USE_PNETCDF_VARN
    PUBLIC USE_PNETCDF_VARN_ON_READ)
target_link_libraries (piof
    $ENV{PNETCDF}/lib/libpnetcdf.a)
