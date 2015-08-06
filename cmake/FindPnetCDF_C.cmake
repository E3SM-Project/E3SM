# - Try to find PnetCDF C Libraries
#
# This can be controlled by setting the PnetCDF_DIR CMake variable, or the
# PNETCDF environment variable.  Alternately, the following CMake variables
# can be set:
#
#   PnetCDF_C_INCLUDE_DIR  (PATH)
#   PnetCDF_C_LIBRARY      (FILE)
#
# Once done, this will define:
#
#   PnetCDF_C_FOUND        (BOOL) - system has NetCDF
#   PnetCDF_C_INCLUDE_DIRS (LIST) - the NetCDF include directories
#   PnetCDF_C_LIBRARIES    (LIST) - link these to use NetCDF
#   PnetCDF_C_DEFINITIONS  (LIST) - preprocessor macros to use with NetCDF
#   PnetCDF_C_OPTIONS      (LIST) - compiler options to use NetCDF

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

# Search for include file
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
if (PnetCDF_C_DIR)
    list (APPEND PnetCDF_C_LIBRARY_HINTS ${PnetCDF_C_DIR}/lib)
endif ()
if (DEFINED ENV{PNETCDF})
    list (APPEND PnetCDF_C_LIBRARY_HINTS $ENV{PNETCDF}/lib)
endif ()

# Search for library file
if (BUILD_SHARED_LIBS)
    find_library (PnetCDF_C_LIBRARY
                  NAMES pnetcdf
                  HINTS ${PnetCDF_C_LIBRARY_HINTS})
else ()
    find_library (PnetCDF_C_LIBRARY
                  NAMES libpnetcdf.a
                  HINTS ${PnetCDF_C_LIBRARY_HINTS})
endif ()

# Unset include search variables
unset (PnetCDF_C_LIBRARY_HINTS)

# Set return variables
set (PnetCDF_C_INCLUDE_DIRS ${PnetCDF_C_INCLUDE_DIR} )
set (PnetCDF_C_LIBRARIES ${PnetCDF_C_LIBRARY} )
set (PnetCDF_C_DEFINITIONS)
set (PnetCDF_C_OPTIONS)

# If static, look for dependencies
if (NOT BUILD_SHARED_LIBS)

    # Dependency find_package arguments
    set (find_args)
    if (PnetCDF_C_FIND_REQUIRED)
        list (APPEND find_args REQUIRED)
    endif ()
    if (PnetCDF_C_FIND_QUIETLY)
        list (APPEND find_args QUIET)
    endif ()

    # DEPENDENCY: MPI
    find_package (MPI ${find_args})
    if (MPI_C_FOUND)
        list (APPEND PnetCDF_C_INCLUDE_DIRS ${MPI_C_INCLUDE_PATH})
        list (APPEND PnetCDF_C_LIBRARIES ${MPI_C_LIBRARIES})
        list (APPEND PnetCDF_C_OPTIONS ${MPI_C_COMPILE_FLAGS})
    endif ()

endif ()

# Check library for varn functions
if (PnetCDF_C_LIBRARY)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_C_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_C_VARN)
    if (PnetCDF_C_VARN)
      list (APPEND PnetCDF_C_DEFINITIONS -DUSE_PNETCDF_VARN)
      list (APPEND PnetCDF_C_DEFINITIONS -DUSE_PNETCDF_VARN_ON_READ)
    endif()  
endif ()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and 
# set PnetCDF_C_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args (PnetCDF_C DEFAULT_MSG
                                   PnetCDF_C_LIBRARY PnetCDF_C_INCLUDE_DIR)
mark_as_advanced (PnetCDF_C_INCLUDE_DIR PnetCDF_C_LIBRARY)
