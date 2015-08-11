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

set (HDF5_VALID_COMPONENTS C HL Fortran Fortran_HL)

if (NOT HDF5_FIND_COMPONENTS)
    set (HDF5_FIND_COMPONENTS C)
endif ()

set (HDF5_FIND_VALID_COMPONENTS)
foreach (comp IN LISTS HDF5_FIND_COMPONENTS)
    if (";${HDF5_VALID_COMPONENTS};" MATCHES ";${comp};")
        list (APPEND HDF5_FIND_VALID_COMPONENTS ${comp})
    endif ()
endforeach ()

set (HDF5_C_INCLUDE_NAMES hdf5.h)
set (HDF5_HL_INCLUDE_NAMES hdf5_hl.h)
set (HDF5_Fortran_INCLUDE_NAMES hdf5.mod)
set (HDF5_Fortran_HL_INCLUDE_NAMES hdf5.mod)

set (HDF5_C_LIBRARY_NAMES hdf5)
set (HDF5_HL_LIBRARY_NAMES hdf5_hl)
set (HDF5_Fortran_LIBRARY_NAMES hdf5_fortran)
set (HDF5_Fortran_HL_LIBRARY_NAMES hdf5hl_fortran)

foreach (comp IN LISTS HDF5_FIND_VALID_COMPONENTS)

    if (NOT HDF5_${comp}_FOUND)

        if (MPI_${comp}_FOUND)
            set (HDF5_${comp}_INCLUDE_HINTS ${MPI_${comp}_INCLUDE_PATH})
            set (HDF5_${comp}_LIBRARY_HINTS)
            foreach (lib IN LISTS ${MPI_${comp}_LIBRARIES})
                get_filename_component (libdir ${lib} PATH)
                list (APPEND HDF5_${comp}_LIBRARY_HINTS ${libdir})
                unset (libdir)
            endforeach ()
        endif ()
    
        find_package_component(HDF5 COMPONENT ${comp}
                               INCLUDE_NAMES ${HDF5_${comp}_INCLUDE_NAMES}
                               INCLUDE_HINTS ${HDF5_${comp}_INCLUDE_HINTS}
                               LIBRARY_NAMES ${HDF5_${comp}_LIBRARY_NAMES}
                               LIBRARY_HINTS ${HDF5_${comp}_LIBRARY_HINTS})
        
        # Handle Dependencies, if static
        if (NOT HDF5_${comp}_IS_SHARED)
        
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
    
        endif ()
        
    endif ()
    
endforeach ()
