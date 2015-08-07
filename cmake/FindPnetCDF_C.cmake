# - Try to find PnetCDF C Libraries
#
# This can be controlled by setting the PnetCDF_DIR CMake variable, or the
# PNETCDF environment variable.
#
# Once done, this will define:
#
#   PnetCDF_C_IS_SHARED    (BOOL) - Whether library is shared/dynamic
#   PnetCDF_C_FOUND        (BOOL) - system has PnetCDF
#   PnetCDF_C_INCLUDE_DIR  (PATH) - Location of the PnetCDF C header
#   PnetCDF_C_INCLUDE_DIRS (LIST) - the NetCDF include directories
#   PnetCDF_C_LIBRARY      (FILE) - Full path to PnetCDF C library file
#   PnetCDF_C_LIBRARIES    (LIST) - link these to use PnetCDF
#   PnetCDF_C_DEFINITIONS  (LIST) - preprocessor macros to use with PnetCDF
#   PnetCDF_C_OPTIONS      (LIST) - compiler options to use PnetCDF


# Determine include dir search order
set (PnetCDF_C_INCLUDE_HINTS)
if (PnetCDF_C_INCLUDE_DIR)
    list (APPEND PnetCDF_C_INCLUDE_HINTS ${PnetCDF_C_INCLUDE_DIR})
endif ()
if (PnetCDF_DIR)
    list (APPEND PnetCDF_C_INCLUDE_HINTS ${PnetCDF_DIR}/include)
endif ()
if (DEFINED ENV{PNETCDF})
    list (APPEND PnetCDF_C_INCLUDE_HINTS $ENV{PNETCDF}/include)
endif ()

# Search for C include file
find_path (PnetCDF_C_INCLUDE_DIR
           NAMES pnetcdf.h
           HINTS ${PnetCDF_C_INCLUDE_HINTS})

# Unset include search variables
unset (PnetCDF_C_INCLUDE_HINTS)

# Determine library dir search order
set (PnetCDF_C_LIBRARY_HINTS)
if (PnetCDF_C_LIBRARY)
    get_filename_component (pnetcdf_library_path ${PnetCDF_C_LIBRARY} PATH)
    list (APPEND PnetCDF_C_LIBRARY_HINTS ${pnetcdf_library_path})
    unset (pnetcdf_library_path)
endif ()
if (PnetCDF_DIR)
    list (APPEND PnetCDF_C_LIBRARY_HINTS ${PnetCDF_DIR}/lib)
endif ()
if (DEFINED ENV{PNETCDF})
    list (APPEND PnetCDF_C_LIBRARY_HINTS $ENV{PNETCDF}/lib)
endif ()

# Search for library file
set (PnetCDF_C_IS_SHARED FALSE)
include (LibFindLibraryMacros)
find_shared_library (PnetCDF_C_LIBRARY
                     NAMES pnetcdf
                     HINTS ${PnetCDF_C_LIBRARY_HINTS}
                     NO_DEFAULT_PATH)
if (PnetCDF_C_LIBRARY)
    set (PnetCDF_C_IS_SHARED TRUE)
else ()
    find_static_library (PnetCDF_C_LIBRARY
                         NAMES pnetcdf
                         HINTS ${PnetCDF_C_LIBRARY_HINTS}
                         NO_DEFAULT_PATH)
endif ()

# Unset include search variables
unset (PnetCDF_C_LIBRARY_HINTS)

# Set C-language return variables
set (PnetCDF_C_INCLUDE_DIRS ${PnetCDF_C_INCLUDE_DIR})
set (PnetCDF_C_LIBRARIES ${PnetCDF_C_LIBRARY})
set (PnetCDF_C_DEFINITIONS)
set (PnetCDF_C_OPTIONS)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and 
# set PnetCDF_C_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args (PnetCDF_C DEFAULT_MSG
                                   PnetCDF_C_LIBRARY PnetCDF_C_INCLUDE_DIR)
mark_as_advanced (PnetCDF_C_INCLUDE_DIR PnetCDF_C_LIBRARY)

# HACK For bug in CMake v3.0:
set (PnetCDF_C_FOUND ${PNETCDF_C_FOUND})
