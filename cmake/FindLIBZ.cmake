# - Try to find LIBZ
#
# This can be controlled by setting the LIBZ_DIR (or, equivalently, the 
# LIBZ environment variable).
#
# Once done, this will define:
#
#   LIBZ_FOUND        (BOOL) - system has LIBZ
#   LIBZ_IS_SHARED    (BOOL) - whether library is shared/dynamic
#   LIBZ_INCLUDE_DIR  (PATH) - Location of the C header file
#   LIBZ_INCLUDE_DIRS (LIST) - the LIBZ include directories
#   LIBZ_LIBRARY      (FILE) - Path to the C library file
#   LIBZ_LIBRARIES    (LIST) - link these to use LIBZ
#
include (LibFind)

# Define LIBZ package
define_package_component (LIBZ
                          INCLUDE_NAMES zlib.h
                          LIBRARY_NAMES z)

# SEARCH FOR PACKAGE
if (NOT LIBZ_FOUND)

    # Manually add the MPI include and library dirs to search paths
    if (MPI_C_FOUND)
        set (LIBZ_PATHS ${MPI_C_INCLUDE_PATH})
        foreach (lib IN LISTS MPI_C_LIBRARIES)
            get_filename_component (libdir ${lib} PATH)
            list (APPEND LIBZ_PATHS ${libdir})
            unset (libdir)
        endforeach ()
    endif ()
    
    # Search for the package
    find_package_component(LIBZ
                           PATHS ${LIBZ_PATHS})

endif ()
