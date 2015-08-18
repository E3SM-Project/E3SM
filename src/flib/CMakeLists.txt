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

set (PIO_Fortran_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
set (PIO_Fortran_DEFINITIONS)
set (PIO_Fortran_COMPILE_FLAGS)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    list (APPEND PIO_Fortran_COMPILE_FLAGS -ffree-line-length-none)
endif()

# Look for c_sizeof capability
check_macro (Fortran_CSIZEOF
                NAME TryCSizeOf.f90
                HINTS ${CMAKE_UTILS_DIR}
                COMMENT "whether Fortran compiler supports c_sizeof")
if (NOT Fortran_CSIZEOF)
    list (APPEND PIO_Fortran_DEFINITIONS NO_C_SIZEOF)
endif()

# Look for filesystem hints
if (PIO_FILESYSTEM_HINTS)
    if (PIO_FILESYSTEM_HINTS STREQUAL lustre)
        message (STATUS "PIO will use lustre filesystem hints")
        list (APPEND PIO_Fortran_DEFINITIONS PIO_LUSTRE_HINTS)
    elseif (PIO_FILESYSTEM_HINTS STREQUAL gpfs)
        message (STATUS "PIO will use gpfs filesystem hints")
        list (APPEND PIO_Fortran_DEFINITIONS PIO_GPFS_HINTS)
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
    list (APPEND PIO_Fortran_INCLUDE_DIRS ${GPTL_INCLUDE_DIRS})
    list (APPEND PIO_Fortran_DEFINITIONS TIMING ${GPTL_DEFINITIONS})
    target_link_libraries (piof gptl)
endif ()

#===== MPI =====
find_package (MPI REQUIRED)

# Check MPI I/O capabilities
if (PIO_ENABLE_MPIIO)
    check_macro (MPIIO_DETECTED
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
    list (APPEND PIO_Fortran_DEFINITIONS USEMPIIO)
endif ()

# Check for MPI Fortran module
check_macro (MPIMOD_DETECTED
                NAME TryMPIMod.f90
                HINTS ${CMAKE_UTILS_DIR}
                COMMENT "whether MPI Fortran module is supported")
if (${MPIMOD_DETECTED})
    message (STATUS "MPI Fortran module detected and enabled.")
else ()
    message (STATUS "MPI Fortran module not detected and therefore disabled.")
    list (APPEND PIO_Fortran_DEFINITIONS NO_MPIMOD)
endif ()

#===== NetCDF-Fortran =====
find_package (NetCDF "4.3.3" COMPONENTS Fortran)
if (NetCDF_Fortran_FOUND)
    list (APPEND PIO_Fortran_INCLUDE_DIRS ${NetCDF_Fortran_INCLUDE_DIRS})
    list (APPEND PIO_Fortran_DEFINITIONS _NETCDF)
    if (NetCDF_C_HAS_PARALLEL)
        list (APPEND PIO_Fortran_DEFINITIONS _NETCDF4)
    endif ()
    target_link_libraries (piof ${NetCDF_Fortran_LIBRARIES})
else ()
    list (APPEND PIO_Fortran_DEFINITIONS _NONETCDF)
endif ()

#===== PnetCDF =====
find_package (PnetCDF "1.6" COMPONENTS Fortran)
if (PnetCDF_Fortran_FOUND)
    list (APPEND PIO_Fortran_INCLUDE_DIRS ${PnetCDF_Fortran_INCLUDE_DIRS})
    list (APPEND PIO_Fortran_DEFINITIONS _PNETCDF)
    target_link_libraries (piof ${PnetCDF_Fortran_LIBRARIES})

    # Check library for varn functions
    include(CheckFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_Fortran_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_VARN)
    if (PnetCDF_VARN)
        list (APPEND PIO_Fortran_DEFINITIONS USE_PNETCDF_VARN USE_PNETCDF_VARN_ON_READ)
    endif()  
else ()
    list (APPEND PIO_Fortran_DEFINITIONS _NOPNETCDF)
endif ()

#===== Add EXTRAs =====
list (APPEND PIO_Fortran_INCLUDE_DIRS ${PIO_Fortran_EXTRA_INCLUDE_DIRS})
list (APPEND PIO_Fortran_DEFINITIONS ${PIO_Fortran_EXTRA_DEFINITIONS})
list (APPEND PIO_Fortran_COMPILE_FLAGS ${PIO_Fortran_EXTRA_COMPILE_FLAGS})
target_link_libraries (piof ${PIO_Fortran_EXTRA_LIBRARIES})

#===== Check for necessities =====
if (NOT PnetCDF_Fortran_FOUND AND NOT NetCDF_Fortran_FOUND)
    message (FATAL_ERROR "Must have PnetCDF and/or NetCDF Fortran libraries")
endif ()

#==============================================================================
#  SET TARGET PROPERTIES
#==============================================================================

set_target_properties (piof PROPERTIES
    INCLUDE_DIRECTORIES "${PIO_Fortran_INCLUDE_DIRS}"
    COMPILE_DEFINITIONS "${PIO_Fortran_DEFINITIONS}"
    COMPILE_FLAGS "${PIO_Fortran_COMPILE_FLAGS}")

# Cache for out-of-scope use
set (PIO_Fortran_INCLUDE_DIRS "${PIO_Fortran_INCLUDE_DIRS}"
     CACHE STRING "PIO Fortran include directories")
set (PIO_Fortran_DEFINITIONS "${PIO_Fortran_DEFINITIONS}"
     CACHE STRING "PIO Fortran compiler definitions")
set (PIO_Fortran_COMPILE_FLAGS "${PIO_Fortran_COMPILE_FLAGS}"
     CACHE STRING "PIO Fortran compiler flags")
     