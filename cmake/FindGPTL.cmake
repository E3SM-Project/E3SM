# - Try to find GPTL
#
# This can be controlled by setting the GPTL_DIR (or, equivalently, the 
# GPTL environment variable), or GPTL_<lang>_DIR CMake variables, where
# <lang> is the COMPONENT language one needs.
#
# Once done, this will define:
#
#   GPTL_<lang>_FOUND        (BOOL) - system has GPTL
#   GPTL_<lang>_IS_SHARED    (BOOL) - whether library is shared/dynamic
#   GPTL_<lang>_INCLUDE_DIR  (PATH) - Location of the C header file
#   GPTL_<lang>_INCLUDE_DIRS (LIST) - the GPTL include directories
#   GPTL_<lang>_LIBRARY      (FILE) - Path to the C library file
#   GPTL_<lang>_LIBRARIES    (LIST) - link these to use GPTL
#   GPTL_<lang>_DEFINITIONS  (LIST) - preprocessor macros to use with GPTL
#   GPTL_<lang>_OPTIONS      (LIST) - compiler options to use GPTL
#
# The available COMPONENTS are: C Fortran Perfmod
# If no components are specified, it assumes only C
include (LibFind)

# Define GPTL C Component
define_package_component (GPTL DEFAULT
                          COMPONENT C
                          INCLUDE_NAMES gptl.h
                          LIBRARY_NAMES gptl)

# Define GPTL Fortran Component
define_package_component (GPTL
                          COMPONENT Fortran
                          INCLUDE_NAMES gptl.mod
                          LIBRARY_NAMES gptl)

# Define GPTL Perfmod Component
define_package_component (GPTL
                          COMPONENT Perfmod
                          INCLUDE_NAMES perf_mod.mod
                          LIBRARY_NAMES gptl)

# Search for list of valid components requested
find_valid_components (GPTL)

#==============================================================================
# SEARCH FOR VALIDATED COMPONENTS
foreach (GPTL_comp IN LISTS GPTL_FIND_VALID_COMPONENTS)

    # If not found already, search...
    if (NOT GPTL_${GPTL_comp}_FOUND)

        # Manually add the MPI include and library dirs to search paths
        if (GPTL_comp STREQUAL C OR GPTL_comp STREQUAL HL)
            if (MPI_C_FOUND)
                set (GPTL_${GPTL_comp}_INCLUDE_HINTS ${MPI_C_INCLUDE_PATH})
                set (GPTL_${GPTL_comp}_LIBRARY_HINTS)
                foreach (lib IN LISTS MPI_C_LIBRARIES)
                    get_filename_component (libdir ${lib} PATH)
                    list (APPEND GPTL_${GPTL_comp}_LIBRARY_HINTS ${libdir})
                    unset (libdir)
                endforeach ()
            endif ()
        else ()
            if (MPI_Fortran_FOUND)
                set (GPTL_${GPTL_comp}_INCLUDE_HINTS ${MPI_Fortran_INCLUDE_PATH})
                set (GPTL_${GPTL_comp}_LIBRARY_HINTS)
                foreach (lib IN LISTS MPI_Fortran_LIBRARIES)
                    get_filename_component (libdir ${lib} PATH)
                    list (APPEND GPTL_${GPTL_comp}_LIBRARY_HINTS ${libdir})
                    unset (libdir)
                endforeach ()
            endif ()
        endif ()

        # Search for the package component
        find_package_component(GPTL COMPONENT ${GPTL_comp}
                               INCLUDE_HINTS ${GPTL_${GPTL_comp}_INCLUDE_HINTS}
                               LIBRARY_HINTS ${GPTL_${GPTL_comp}_LIBRARY_HINTS})

    endif ()
    
endforeach ()
