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

# SEARCH FOR VALIDATED COMPONENTS
foreach (comp IN LISTS HDF5_FIND_VALID_COMPONENTS)

    # If not found already, search...
    if (NOT HDF5_${comp}_FOUND)

        # Manually add the MPI include and library dirs to search paths
        if (MPI_${comp}_FOUND)
            set (HDF5_${comp}_INCLUDE_HINTS ${MPI_${comp}_INCLUDE_PATH})
            set (HDF5_${comp}_LIBRARY_HINTS)
            foreach (lib IN LISTS MPI_${comp}_LIBRARIES)
                get_filename_component (libdir ${lib} PATH)
                list (APPEND HDF5_${comp}_LIBRARY_HINTS ${libdir})
                unset (libdir)
            endforeach ()
        endif ()
    
        # Search for the package component
        find_package_component(HDF5 COMPONENT ${comp}
                               INCLUDE_NAMES ${HDF5_${comp}_INCLUDE_NAMES}
                               INCLUDE_HINTS ${HDF5_${comp}_INCLUDE_HINTS}
                               LIBRARY_NAMES ${HDF5_${comp}_LIBRARY_NAMES}
                               LIBRARY_HINTS ${HDF5_${comp}_LIBRARY_HINTS})
        
    endif ()
    
endforeach ()

# SEARCH FOR DEPENDENCIES (only if SHARED libraries were found)
foreach (comp IN LISTS HDF5_FIND_VALID_COMPONENTS)

    # If the component was found, and it is a static library...
    if (HDF5_${comp}_FOUND AND NOT HDF5_${comp}_IS_SHARED)
        
        # Search only if dependencies for this component were not already found
        if (NOT HDF5_${comp}_DEPENDENCIES_SEARCHED)

            # COMPONENT: C
            if (comp STREQUAL C)

                # DEPENDENCY: LIBZ
                find_package (LIBZ)
                if (LIBZ_FOUND)
                    list (APPEND HDF5_C_INCLUDE_DIRS ${LIBZ_INCLUDE_DIRS})
                    list (APPEND HDF5_C_LIBRARIES ${LIBZ_LIBRARIES})
                endif ()

            else ()

                # DEPENDENCY: HDF5
                find_package (HDF5 COMPONENTS C)
                if (HDF5_C_FOUND)
                    list (APPEND HDF5_${comp}_INCLUDE_DIRS ${HDF5_C_INCLUDE_DIRS})
                    list (APPEND HDF5_${comp}_LIBRARIES ${HDF5_C_LIBRARIES})
                endif ()
    
            endif ()
            
            set (HDF5_${comp}_DEPENDENCIES_SEARCHED TRUE)
    
        endif ()
        
    endif ()
    
endforeach ()
