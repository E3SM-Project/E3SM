cmake_minimum_required (VERSION 2.8.10)
project (PIOF Fortran)
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

add_library (piof ${PIO_Fortran_SRCS} ${PIO_GenF90_SRCS})
if (NOT PIO_ENABLE_FORTRAN)
    set_target_properties(piof PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif ()

include_directories (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set_target_properties(piof PROPERTIES COMPILE_FLAGS -ffree-line-length-none)
endif()

# Look for c_sizeof capability
check_property (Fortran_CSIZEOF
                NAME TryCSizeOf.f90
                HINTS ${CMAKE_UTILS_DIR}
                COMMENT "whether Fortran compiler supports c_sizeof")
if (Fortran_CSIZEOF)
    add_definitions (-DNO_C_SIZEOF)
endif()

# Look for filesystem hints
if (PIO_FILESYSTEM_HINTS)
    if (PIO_FILESYSTEM_HINTS STREQUAL lustre)
        message (STATUS "PIO will use lustre filesystem hints")
        add_definitions (-DPIO_LUSTRE_HINTS)
    elseif (PIO_FILESYSTEM_HINTS STREQUAL gpfs)
        message (STATUS "PIO will use gpfs filesystem hints")
        add_definitions (-DPIO_GPFS_HINTS)
    else ()
        message (WARNING "${PIO_FILESYSTEM_HINTS} not valid option for PIO_FILESYSTEM_HINTS; use gpfs or lustre.")
    endif ()
endif ()

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================

#===== genf90 =====
if (DEFINED GENF90_DIR)
    add_custom_target(genf90
        DEPENDS ${GENF90_DIR}/genf90.pl)
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
    set (GENF90_DIR ${SOURCE_DIR})
    unset (SOURCE_DIR)
endif ()
add_dependencies (piof genf90)

#===== Fortran Source Generation with GenF90 =====
foreach (SRC_FILE IN LISTS PIO_GenF90_SRCS)
    add_custom_command (OUTPUT ${SRC_FILE}
        COMMAND ${GENF90_DIR}/genf90.pl
                ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in > ${SRC_FILE}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in genf90)
endforeach ()

#===== GPTL =====
if (PIO_ENABLE_TIMING)
    add_definitions (-DTIMING)
    include_directories(${GPTL_INCLUDE_DIRS})
    target_link_libraries (piof gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)

# Check MPI I/O capabilities
if (PIO_ENABLE_MPIIO)
    check_property (MPIIO_DETECTED
                    NAME TryMPIIO.f90
                    HINTS ${CMAKE_UTILS_DIR}
                    COMMENT "whether MPIIO is supported")
    if (${MPIIO_DETECTED})
        message (STATUS "MPIIO detected and enabled.")
    else ()
        message (STATUS "MPIIO not detected and therefore disabled.")
        set (PIO_ENABLE_MPIIO FALSE)
    endif ()
endif ()
if (${PIO_ENABLE_MPIIO}) 
    add_definitions (-DUSEMPIIO)
endif ()

# Check for MPI Fortran module
check_property (MPIMOD_DETECTED
                NAME TryMPIMod.f90
                HINTS ${CMAKE_UTILS_DIR}
                COMMENT "whether MPI Fortran module is supported")
if (${MPIMOD_DETECTED})
    message (STATUS "MPI Fortran module detected and enabled.")
else ()
    message (STATUS "MPI Fortran module not detected and therefore disabled.")
    add_definitions (-DNO_MPIMOD)
endif ()

#===== NetCDF-Fortran =====
find_package (NetCDF "4.3.3" COMPONENTS Fortran)
if (NetCDF_Fortran_FOUND)
    include_directories (${NetCDF_Fortran_INCLUDE_DIRS})
    add_definitions (-D_NETCDF)
    target_link_libraries (piof ${NetCDF_Fortran_LIBRARIES})
    if (EXISTS ${NetCDF_Fortran_INCLUDE_DIR}/netcdf_par.h)
        add_definitions (-D_NETCDF4)
    endif ()
else ()
    add_definitions (-D_NONETCDF)
endif ()

#===== PnetCDF =====
find_package (PnetCDF "1.6" COMPONENTS Fortran)
if (PnetCDF_Fortran_FOUND)
    include_directories (${PnetCDF_Fortran_INCLUDE_DIRS})
    add_definitions (-D_PNETCDF)
    target_link_libraries (piof ${PnetCDF_Fortran_LIBRARIES})

    # Check library for varn functions
    include(CheckFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_Fortran_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
        add_definitions(-DUSE_PNETCDF_VARN -DUSE_PNETCDF_VARN_ON_READ)
    endif()  
else ()
    add_definitions (-D_NOPNETCDF)
endif ()

#===== Add EXTRAs =====
include_directories (${PIO_Fortran_EXTRA_INCLUDE_DIRS})
target_link_libraries (pioc ${PIO_Fortran_EXTRA_LIBRARIES})
add_definitions (${PIO_Fortran_EXTRA_DEFINITIONS})
if (PIO_Fortran_EXTRA_LINK_FLAGS)
    set_target_properties(pioc PROPERTIES 
        LINK_FLAGS ${PIO_Fortran_EXTRA_LINK_FLAGS})
endif ()

#===== Check for necessities =====
if (NOT PnetCDF_Fortran_FOUND AND NOT NetCDF_Fortran_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF Fortran libraries")
endif ()

#==============================================================================
#  DEFINE PARAMETERS FOR TARGETS DEPENDENT ON THIS
#==============================================================================

get_target_property(PIO_Fortran_INCLUDE_DIRS piof INCLUDE_DIRECTORIES)
set (PIO_Fortran_INCLUDE_DIRS ${PIO_Fortran_INCLUDE_DIRS}
     CACHE STRING "PIO Fortran Include Directory")
