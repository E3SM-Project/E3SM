cmake_minimum_required (VERSION 2.8.12)
project (PIOF Fortran)
include (CheckFunctionExists)
include (ExternalProject)
include (LibCheck)

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
                     
set (PIO_Fortran_MODS ${CMAKE_CURRENT_BINARY_DIR}/pio.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/pio_nf.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/pio_types.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/piolib_mod.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/pionfget_mod.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/pio_kinds.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/pio_support.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/piodarray.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/pionfatt_mod.mod
                      ${CMAKE_CURRENT_BINARY_DIR}/pionfput_mod.mod)

add_library (piof ${PIO_Fortran_SRCS} ${PIO_GenF90_SRCS})
if (NOT PIO_ENABLE_FORTRAN)
    set_target_properties(piof PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif ()

target_include_directories (piof
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options (piof
        PRIVATE -ffree-line-length-none)
endif()

# Look for c_sizeof capability
check_macro (Fortran_CSIZEOF
             NAME TryCSizeOf.f90
             HINTS ${CMAKE_MODULE_PATH}
             COMMENT "whether Fortran compiler supports c_sizeof")
if (NOT Fortran_CSIZEOF)
    target_compile_definitions (piof
        PUBLIC NO_C_SIZEOF)
endif()

# Look for filesystem hints
if (PIO_FILESYSTEM_HINTS)
    if (PIO_FILESYSTEM_HINTS STREQUAL lustre)
        message (STATUS "PIO will use lustre filesystem hints")
        target_compile_definitions (piof
            PUBLIC PIO_LUSTRE_HINTS)
    elseif (PIO_FILESYSTEM_HINTS STREQUAL gpfs)
        message (STATUS "PIO will use gpfs filesystem hints")
        target_compile_definitions (piof
            PUBLIC PIO_GPFS_HINTS)
    else ()
        message (WARNING "${PIO_FILESYSTEM_HINTS} not valid option for PIO_FILESYSTEM_HINTS; use gpfs or lustre.")
    endif ()
endif ()

#==============================================================================
#  DEFINE THE INSTALL
#==============================================================================

# Library
install (TARGETS piof DESTINATION lib)

# Fortran Modules
install (FILES ${PIO_Fortran_MODS} DESTINATION include)

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================

#===== pioc =====
target_link_libraries(piof
    PUBLIC pioc)
    
#===== genf90 =====
if (DEFINED GENF90_PATH)
    add_custom_target(genf90
        DEPENDS ${GENF90_PATH}/genf90.pl)
else ()
    ExternalProject_Add (genf90
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/genf90
        GIT_REPOSITORY https://github.com/PARALLELIO/genf90
        GIT_TAG genf90_140121
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND "")
    ExternalProject_Get_Property (genf90 SOURCE_DIR)
    set (GENF90_PATH ${SOURCE_DIR})
    unset (SOURCE_DIR)
endif ()
add_dependencies (piof genf90)

#===== Fortran Source Generation with GenF90 =====
foreach (SRC_FILE IN LISTS PIO_GenF90_SRCS)
    add_custom_command (OUTPUT ${SRC_FILE}
        COMMAND ${GENF90_PATH}/genf90.pl
                ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in > ${SRC_FILE}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in genf90)
endforeach ()

#===== MPI =====
find_package (MPI REQUIRED)

# Check MPI I/O capabilities
if (PIO_USE_MPIIO)
    check_macro (MPIIO_DETECTED
                 NAME TryMPIIO.f90
                 HINTS ${CMAKE_MODULE_PATH}
                 COMMENT "whether MPIIO is supported")
    if (${MPIIO_DETECTED})
        message (STATUS "MPIIO detected and enabled.")
    else ()
        message (STATUS "MPIIO not detected and therefore disabled.")
        set (PIO_USE_MPIIO FALSE)
    endif ()
endif ()
if (${PIO_USE_MPIIO}) 
    target_compile_definitions (piof
        PUBLIC USEMPIIO)
endif ()

# Check for MPI Fortran module
check_macro (MPIMOD_DETECTED
             NAME TryMPIMod.f90
             HINTS ${CMAKE_MODULE_PATH}
             COMMENT "whether MPI Fortran module is supported")
if (${MPIMOD_DETECTED})
    message (STATUS "MPI Fortran module detected and enabled.")
else ()
    message (STATUS "MPI Fortran module not detected and therefore disabled.")
    target_compile_definitions (piof
        PUBLIC NO_MPIMOD)
endif ()

#===== GPTL =====
if (PIO_ENABLE_TIMING)
    find_package (GPTL COMPONENTS Fortran_Perf QUIET)
    if (GPTL_Fortran_Perf_FOUND)
        message (STATUS "Found GPTL Fortran Perf: ${GPTL_Fortran_Perf_LIBRARIES}")
        target_include_directories (piof 
            PUBLIC ${GPTL_Fortran_INCLUDE_DIRS})
        target_link_libraries (piof
            PUBLIC ${GPTL_Fortran_LIBRARIES})
    else ()
        message (STATUS "Using internal GPTL Fortran library for timing")
        target_link_libraries (piof
            PUBLIC gptl)
    endif ()
endif ()

#===== NetCDF-Fortran =====
find_package (NetCDF "4.3.3" COMPONENTS Fortran)
if (NetCDF_Fortran_FOUND)
    target_include_directories (piof 
        PUBLIC ${NetCDF_Fortran_INCLUDE_DIRS})
    target_compile_definitions (piof 
        PUBLIC _NETCDF)
    target_link_libraries (piof
        PUBLIC ${NetCDF_Fortran_LIBRARIES})
    if (EXISTS ${NetCDF_Fortran_INCLUDE_DIR}/netcdf_par.h)
        target_compile_definitions (piof
            PUBLIC _NETCDF4)
    endif ()
else ()
    target_compile_definitions (piof 
        PUBLIC _NONETCDF)
endif ()

#===== PnetCDF =====
if (WITH_PNETCDF)
    find_package (PnetCDF "1.6" COMPONENTS Fortran REQUIRED)
endif ()
if (PnetCDF_Fortran_FOUND)
    target_include_directories (piof 
        PUBLIC ${PnetCDF_Fortran_INCLUDE_DIRS})
    target_compile_definitions (piof 
        PUBLIC _PNETCDF)
    target_link_libraries (piof
        PUBLIC ${PnetCDF_Fortran_LIBRARIES})

    # Check library for varn functions
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_Fortran_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
        target_compile_definitions(piof
            PUBLIC USE_PNETCDF_VARN
            PUBLIC USE_PNETCDF_VARN_ON_READ)
    endif()  
else ()
    target_compile_definitions (piof 
        PUBLIC _NOPNETCDF)
endif ()

#===== Add EXTRAs =====
target_include_directories (piof
    PUBLIC ${PIO_Fortran_EXTRA_INCLUDE_DIRS})
target_link_libraries (piof
    PUBLIC ${PIO_Fortran_EXTRA_LIBRARIES})
target_compile_options (piof
    PRIVATE ${PIO_Fortran_EXTRA_COMPILE_OPTIONS})
target_compile_definitions (piof
    PUBLIC ${PIO_Fortran_EXTRA_COMPILE_DEFINITIONS})
if (PIO_Fortran_EXTRA_LINK_FLAGS)
    set_target_properties(piof PROPERTIES 
        LINK_FLAGS ${PIO_Fortran_EXTRA_LINK_FLAGS})
endif ()

#===== Check for necessities =====
if (NOT PnetCDF_Fortran_FOUND AND NOT NetCDF_Fortran_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF Fortran libraries")
endif ()
     