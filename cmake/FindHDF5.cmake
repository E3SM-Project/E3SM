# - Try to find HDF5
#
# This can be controlled by setting the HDF5_DIR (or, equivalently, the 
# HDF5 environment variable), or HDF5_<lang>_DIR CMake variables, where
# <lang> is the COMPONENT language one needs.
#
# Once done, this will define:
#
#   HDF5_<lang>_FOUND        (BOOL) - system has HDF5
#   HDF5_<lang>_IS_SHARED    (BOOL) - whether library is shared/dynamic
#   HDF5_<lang>_INCLUDE_DIR  (PATH) - Location of the C header file
#   HDF5_<lang>_INCLUDE_DIRS (LIST) - the HDF5 include directories
#   HDF5_<lang>_LIBRARY      (FILE) - Path to the C library file
#   HDF5_<lang>_LIBRARIES    (LIST) - link these to use HDF5
#   HDF5_<lang>_DEFINITIONS  (LIST) - preprocessor macros to use with HDF5
#   HDF5_<lang>_OPTIONS      (LIST) - compiler options to use HDF5
#
# The available COMPONENTS are: C HL Fortran Fortran_HL
# If no components are specified, it assumes only C
include (LibFindLibraryMacros)
include (CheckHDF5)

# Define HDF5 C Component
define_package_component (HDF5 DEFAULT
                          COMPONENT C
                          INCLUDE_NAMES hdf5.h
                          LIBRARY_NAMES hdf5)

# Define HDF5 HL Component
define_package_component (HDF5
                          COMPONENT HL
                          INCLUDE_NAMES hdf5_hl.h
                          LIBRARY_NAMES hdf5_hl)

# Define HDF5 Fortran Component
define_package_component (HDF5
                          COMPONENT Fortran
                          INCLUDE_NAMES hdf5.mod
                          LIBRARY_NAMES hdf5_fortran)

# Define HDF5 Fortran_HL Component
define_package_component (HDF5
                          COMPONENT Fortran_HL
                          INCLUDE_NAMES hdf5.mod
                          LIBRARY_NAMES hdf5hl_fortran)

# Search for list of valid components requested
find_valid_components (HDF5)

#==============================================================================
# SEARCH FOR VALIDATED COMPONENTS
foreach (HDF5_comp IN LISTS HDF5_FIND_VALID_COMPONENTS)

    # If not found already, search...
    if (NOT HDF5_${HDF5_comp}_FOUND)

        # Manually add the MPI include and library dirs to search paths
        if (HDF5_comp STREQUAL C OR HDF5_comp STREQUAL HL)
            if (MPI_C_FOUND)
                set (HDF5_${HDF5_comp}_INCLUDE_HINTS ${MPI_C_INCLUDE_PATH})
                set (HDF5_${HDF5_comp}_LIBRARY_HINTS)
                foreach (lib IN LISTS MPI_C_LIBRARIES)
                    get_filename_component (libdir ${lib} PATH)
                    list (APPEND HDF5_${HDF5_comp}_LIBRARY_HINTS ${libdir})
                    unset (libdir)
                endforeach ()
            endif ()
        else ()
            if (MPI_Fortran_FOUND)
                set (HDF5_${HDF5_comp}_INCLUDE_HINTS ${MPI_Fortran_INCLUDE_PATH})
                set (HDF5_${HDF5_comp}_LIBRARY_HINTS)
                foreach (lib IN LISTS MPI_Fortran_LIBRARIES)
                    get_filename_component (libdir ${lib} PATH)
                    list (APPEND HDF5_${HDF5_comp}_LIBRARY_HINTS ${libdir})
                    unset (libdir)
                endforeach ()
            endif ()
        endif ()

        # Search for the package component
        find_package_component(HDF5 COMPONENT ${HDF5_comp}
                               INCLUDE_HINTS ${HDF5_${HDF5_comp}_INCLUDE_HINTS}
                               LIBRARY_HINTS ${HDF5_${HDF5_comp}_LIBRARY_HINTS})

        # Dependencies
        if (HDF5_comp STREQUAL C AND NOT HDF5_C_IS_SHARED)
        
            # DEPENDENCY: LIBZ
            find_package (LIBZ)
            if (LIBZ_FOUND)
                list (APPEND HDF5_C_INCLUDE_DIRS ${LIBZ_INCLUDE_DIRS})
                list (APPEND HDF5_C_LIBRARIES ${LIBZ_LIBRARIES})
            endif ()
            
        elseif (NOT HDF5_${HDF5_comp}_IS_SHARED)

            # DEPENDENCY: HDF5
            find_package (HDF5 COMPONENTS C)
            if (HDF5_C_FOUND)
                list (APPEND HDF5_${HDF5_comp}_INCLUDE_DIRS ${HDF5_C_INCLUDE_DIRS})
                list (APPEND HDF5_${HDF5_comp}_LIBRARIES ${HDF5_C_LIBRARIES})
            endif ()

        endif ()

    endif ()
    
endforeach ()

#==============================================================================
# CHECKS AND DEPENDENCIES
foreach (HDF5_comp IN LISTS HDF5_FIND_VALID_COMPONENTS)
    if (HDF5_comp STREQUAL C)
        check_HDF5_C ()
    endif ()
endforeach ()
