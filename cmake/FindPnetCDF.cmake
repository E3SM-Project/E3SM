# - Try to find PnetCDF
#
# This can be controlled by setting the PnetCDF_DIR (or, equivalently, the 
# PNETCDF environment variable), or PnetCDF_<lang>_DIR CMake variables, where
# <lang> is the COMPONENT language one needs.
#
# Once done, this will define:
#
#   PnetCDF_<lang>_FOUND        (BOOL) - system has PnetCDF
#   PnetCDF_<lang>_IS_SHARED    (BOOL) - whether library is shared/dynamic
#   PnetCDF_<lang>_INCLUDE_DIR  (PATH) - Location of the C header file
#   PnetCDF_<lang>_INCLUDE_DIRS (LIST) - the PnetCDF include directories
#   PnetCDF_<lang>_LIBRARY      (FILE) - Path to the C library file
#   PnetCDF_<lang>_LIBRARIES    (LIST) - link these to use PnetCDF
#   PnetCDF_<lang>_DEFINITIONS  (LIST) - preprocessor macros to use with PnetCDF
#   PnetCDF_<lang>_OPTIONS      (LIST) - compiler options to use PnetCDF
#
# The available COMPONENTS are: C, Fortran
# If no components are specified, it assumes only C
include (LibFindLibraryMacros)

# Define PnetCDF C Component
define_package_component (PnetCDF DEFAULT
                          COMPONENT C
                          INCLUDE_NAMES pnetcdf.h
                          LIBRARY_NAMES pnetcdf)

# Define PnetCDF Fortran Component
define_package_component (PnetCDF
                          COMPONENT Fortran
                          INCLUDE_NAMES pnetcdf.mod pnetcdf.inc
                          LIBRARY_NAMES pnetcdf)

# Search for list of valid components requested
find_valid_components (PnetCDF)

# SEARCH FOR VALIDATED COMPONENTS
foreach (PnetCDF_comp IN LISTS PnetCDF_FIND_VALID_COMPONENTS)

    # If not found already, search...
    if (NOT PnetCDF_${PnetCDF_comp}_FOUND)

        # Manually add the MPI include and library dirs to search paths
        if (MPI_${PnetCDF_comp}_FOUND)
            set (PnetCDF_${PnetCDF_comp}_INCLUDE_HINTS ${MPI_${PnetCDF_comp}_INCLUDE_PATH})
            set (PnetCDF_${PnetCDF_comp}_LIBRARY_HINTS)
            foreach (lib IN LISTS MPI_${PnetCDF_comp}_LIBRARIES)
                get_filename_component (libdir ${lib} PATH)
                list (APPEND PnetCDF_${PnetCDF_comp}_LIBRARY_HINTS ${libdir})
                unset (libdir)
            endforeach ()
        endif ()
        
        # Search for the package component
        find_package_component(PnetCDF COMPONENT ${PnetCDF_comp}
                               INCLUDE_HINTS ${PnetCDF_${PnetCDF_comp}_INCLUDE_HINTS}
                               LIBRARY_HINTS ${PnetCDF_${PnetCDF_comp}_LIBRARY_HINTS})
    
        # Continue only if found
        if (PnetCDF_${PnetCDF_comp}_FOUND)

        #----------------------------------------------------------------------
        # Checks & Dependencies for COMPONENT: C
        if (PnetCDF_comp STREQUAL C AND NOT PnetCDF_C_FINISHED)

            # Get version string
            try_run (PnetCDF_C_VERSION_RUNVAR PnetCDF_C_VERSION_COMPVAR
                     ${CMAKE_CURRENT_BINARY_DIR}/tryPnetCDF_C_VERSION
                     ${CMAKE_SOURCE_DIR}/cmake/TryPnetCDF_VERSION.c
                     COMPILE_DEFINITIONS -I${PnetCDF_C_INCLUDE_DIR}
                     COMPILE_OUTPUT_VARIABLE TryPnetCDF_OUT
                     RUN_OUTPUT_VARIABLE PnetCDF_C_VERSION)
            if (PnetCDF_C_VERSION)
                if (PnetCDF_C_VERSION VERSION_LESS PnetCDF_FIND_VERSION)
                    message (FATAL_ERROR "PnetCDF_C version insufficient")
                else ()
                    message (STATUS "Found PnetCDF_C version ${PnetCDF_C_VERSION}")
                endif ()
            else ()
                message (STATUS "Could not find PnetCDF_C version")
            endif ()

            # Checks and dependencies finished
            set (PnetCDF_C_FINISHED TRUE
                 CACHE BOOL "PnetCDF C fully found")

        #----------------------------------------------------------------------
        # Checks & Dependencies for COMPONENT: Fortran
        elseif (PnetCDF_comp STREQUAL Fortran AND NOT PnetCDF_Fortran_FINISHED)

            # Get version string
            set (COMP_DEFS)
            foreach (incdir IN LISTS PnetCDF_Fortran_INCLUDE_DIRS)
                list (APPEND COMP_DEFS "-I${incdir}")
            endforeach ()
            try_run (PnetCDF_Fortran_VERSION_RUNVAR PnetCDF_Fortran_VERSION_COMPVAR
                     ${CMAKE_CURRENT_BINARY_DIR}/tryPnetCDF_Fortran_VERSION
                     ${CMAKE_SOURCE_DIR}/cmake/TryPnetCDF_VERSION.f90
                     COMPILE_DEFINITIONS ${COMP_DEFS}
                     LINK_LIBRARIES ${PnetCDF_Fortran_LIBRARIES}
                     COMPILE_OUTPUT_VARIABLE TryPnetCDF_OUT
                     RUN_OUTPUT_VARIABLE PnetCDF_Fortran_VERSION)
            if (PnetCDF_Fortran_VERSION)
                string (STRIP ${PnetCDF_Fortran_VERSION} PnetCDF_Fortran_VERSION)
                if (PnetCDF_Fortran_VERSION VERSION_LESS PnetCDF_FIND_VERSION)
                    message (FATAL_ERROR "PnetCDF_Fortran version insufficient")
                else ()
                    message (STATUS "Found PnetCDF_Fortran version ${PnetCDF_Fortran_VERSION}")
                endif ()
            else ()
                message (STATUS "Could not find PnetCDF_Fortran version")
            endif ()

            # Checks and dependencies finished
            set (PnetCDF_Fortran_FINISHED TRUE
                 CACHE BOOL "PnetCDF Fortran fully found")

        endif ()

        endif ()
        
    endif ()
    
endforeach ()
