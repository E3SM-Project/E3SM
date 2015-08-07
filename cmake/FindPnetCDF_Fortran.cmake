# - Try to find PnetCDF Fortran Libraries
#
# This can be controlled by setting the PnetCDF_DIR CMake variable, or the
# PNETCDF environment variable.
#
# Once done, this will define:
#
#   PnetCDF_Fortran_IS_SHARED    (BOOL) - whether the PnetCDF library is shared/dynamic
#
#   PnetCDF_Fortran_FOUND        (BOOL) - system has PnetCDF
#   PnetCDF_Fortran_INCLUDE_DIR  (PATH) - Location of the PnetCDF Fortran header
#   PnetCDF_Fortran_INCLUDE_DIRS (LIST) - the PnetCDF include directories
#   PnetCDF_Fortran_LIBRARY      (FILE) - Full path to PnetCDF Fortran library file
#   PnetCDF_Fortran_LIBRARIES    (LIST) - link these to use PnetCDF
#   PnetCDF_Fortran_DEFINITIONS  (LIST) - preprocessor macros to use with PnetCDF
#   PnetCDF_Fortran_OPTIONS      (LIST) - compiler options to use PnetCDF


# Determine include dir search order
set (PnetCDF_Fortran_INCLUDE_HINTS)
if (PnetCDF_Fortran_INCLUDE_DIR)
    list (APPEND PnetCDF_Fortran_INCLUDE_HINTS ${PnetCDF_Fortran_INCLUDE_DIR})
endif ()
if (PnetCDF_DIR)
    list (APPEND PnetCDF_Fortran_INCLUDE_HINTS ${PnetCDF_DIR}/include)
endif ()
if (DEFINED ENV{PNETCDF})
    list (APPEND PnetCDF_Fortran_INCLUDE_HINTS $ENV{PNETCDF}/include)
endif ()

# Search for Fortran include file
find_path (PnetCDF_Fortran_INCLUDE_DIR
           NAMES pnetcdf.inc pnetcdf.mod
           HINTS ${PnetCDF_Fortran_INCLUDE_HINTS})
           
# Unset include search variables
unset (PnetCDF_Fortran_INCLUDE_HINTS)

# Determine library dir search order
set (PnetCDF_Fortran_LIBRARY_HINTS)
if (PnetCDF_Fortran_LIBRARY)
    get_filename_component (pnetcdf_library_path ${PnetCDF_Fortran_LIBRARY} PATH)
    list (APPEND PnetCDF_Fortran_LIBRARY_HINTS ${pnetcdf_library_path})
    unset (pnetcdf_library_path)
endif ()
if (PnetCDF_DIR)
    list (APPEND PnetCDF_Fortran_LIBRARY_HINTS ${PnetCDF_DIR}/lib)
endif ()
if (DEFINED ENV{PNETCDF})
    list (APPEND PnetCDF_Fortran_LIBRARY_HINTS $ENV{PNETCDF}/lib)
endif ()

# Search for shared and static library files
set (PnetCDF_Fortran_IS_SHARED FALSE)
find_library (PnetCDF_Fortran_LIBRARY
              NAMES pnetcdff
              HINTS ${PnetCDF_Fortran_LIBRARY_HINTS})
if (PnetCDF_Fortran_LIBRARY)
    set (PnetCDF_Fortran_IS_SHARED TRUE)
else ()
    find_library (PnetCDF_Fortran_LIBRARY
                  NAMES libpnetcdff.a
                  HINTS ${PnetCDF_Fortran_LIBRARY_HINTS})
endif ()

# Unset include search variables
unset (PnetCDF_Fortran_LIBRARY_HINTS)

# Set Fortran-language return variables
set (PnetCDF_Fortran_INCLUDE_DIRS ${PnetCDF_Fortran_INCLUDE_DIR})
set (PnetCDF_Fortran_LIBRARIES ${PnetCDF_Fortran_LIBRARY})
set (PnetCDF_Fortran_DEFINITIONS)
set (PnetCDF_Fortran_OPTIONS)

# If static, look for dependencies
if (NOT PnetCDF_Fortran_IS_SHARED)

    # Dependency find_package arguments
    set (find_args)
    if (PnetCDF_Fortran_FIND_REQUIRED)
        list (APPEND find_args REQUIRED)
    endif ()
    if (PnetCDF_Fortran_FIND_QUIETLY)
        list (APPEND find_args QUIET)
    endif ()

    # DEPENDENCY: MPI
    find_package (MPI ${find_args})
    if (MPI_Fortran_FOUND)
        list (APPEND PnetCDF_Fortran_INCLUDE_DIRS ${MPI_Fortran_INCLUDE_PATH})
        list (APPEND PnetCDF_Fortran_LIBRARIES ${MPI_Fortran_LIBRARIES})
        list (APPEND PnetCDF_Fortran_OPTIONS ${MPI_Fortran_COMPILE_FLAGS})
    endif ()

endif ()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and 
# set PnetCDF_Fortran_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args (PnetCDF_Fortran DEFAULT_MSG
                                   PnetCDF_Fortran_LIBRARY PnetCDF_Fortran_INCLUDE_DIR)
mark_as_advanced (PnetCDF_Fortran_INCLUDE_DIR PnetCDF_Fortran_LIBRARY)

# HACK For bug in CMake v3.0:
set (PnetCDF_Fortran_FOUND ${PNETCDF_FORTRAN_FOUND})
