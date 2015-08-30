# - Try to find SZIP
#
# This can be controlled by setting the SZIP_DIR (or, equivalently, the 
# SZIP environment variable).
#
# Once done, this will define:
#
#   SZIP_FOUND        (BOOL) - system has SZIP
#   SZIP_IS_SHARED    (BOOL) - whether library is shared/dynamic
#   SZIP_INCLUDE_DIR  (PATH) - Location of the C header file
#   SZIP_INCLUDE_DIRS (LIST) - the SZIP include directories
#   SZIP_LIBRARY      (FILE) - Path to the C library file
#   SZIP_LIBRARIES    (LIST) - link these to use SZIP
#   SZIP_DEFINITIONS  (LIST) - preprocessor macros to use with SZIP
#   SZIP_OPTIONS      (LIST) - compiler options to use SZIP
#
include (LibFind)

# Define SZIP package
define_package_component (SZIP
                          INCLUDE_NAMES szlib.h
                          LIBRARY_NAMES sz)

# SEARCH FOR PACKAGE
if (NOT SZIP_FOUND)

    # Manually add the MPI include and library dirs to search paths
    if (MPI_C_FOUND)
        set (SZIP_INCLUDE_PATHS ${MPI_C_INCLUDE_PATH})
        set (SZIP_LIBRARY_PATHS)
        foreach (lib IN LISTS MPI_C_LIBRARIES)
            get_filename_component (libdir ${lib} PATH)
            list (APPEND SZIP_LIBRARY_PATHS ${libdir})
            unset (libdir)
        endforeach ()
    endif ()
    
    # Search for the package
    find_package_component(SZIP
                           INCLUDE_PATHS ${SZIP_INCLUDE_PATHS}
                           LIBRARY_PATHS ${SZIP_LIBRARY_PATHS})

endif ()
