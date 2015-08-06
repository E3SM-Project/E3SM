# - Try to find PnetCDF Fortran Libraries
#
# This can be controlled by setting the PnetCDF_DIR CMake variable, or the
# PNETCDF environment variable.  Alternately, the following CMake variables
# can be set:
#
#   PnetCDF_Fortran_INCLUDE_DIR  (PATH)
#   PnetCDF_Fortran_LIBRARY      (FILE)
#
# Once done, this will define:
#
#   PnetCDF_Fortran_FOUND        (BOOL) - system has NetCDF
#   PnetCDF_Fortran_INCLUDE_DIRS (LIST) - the NetCDF include directories
#   PnetCDF_Fortran_LIBRARIES    (LIST) - link these to use NetCDF
#   PnetCDF_Fortran_DEFINITIONS  (LIST) - preprocessor macros to use with NetCDF
#   PnetCDF_Fortran_OPTIONS      (LIST) - compiler options to use NetCDF

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

# Search for include file
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

# Search for library file
if (BUILD_SHARED_LIBS)
    find_library (PnetCDF_Fortran_LIBRARY
                  NAMES pnetcdf
                  HINTS ${PnetCDF_Fortran_LIBRARY_HINTS})
else ()
    find_library (PnetCDF_Fortran_LIBRARY
                  NAMES libpnetcdf.a
                  HINTS ${PnetCDF_Fortran_LIBRARY_HINTS})
endif ()

# Unset include search variables
unset (PnetCDF_Fortran_LIBRARY_HINTS)

# Set return variables
set (PnetCDF_Fortran_INCLUDE_DIRS ${PnetCDF_Fortran_INCLUDE_DIR} )
set (PnetCDF_Fortran_LIBRARIES ${PnetCDF_Fortran_LIBRARY} )
set (PnetCDF_Fortran_DEFINITIONS)
set (PnetCDF_Fortran_OPTIONS)

# If static, look for dependencies
if (NOT BUILD_SHARED_LIBS)

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

# Check library for varn functions
include(CheckFunctionExists)
if (PnetCDF_Fortran_LIBRARY)
    set (CMAKE_REQUIRED_LIBRARIES ${PnetCDF_Fortran_LIBRARY})
    check_function_exists (ncmpi_get_varn PnetCDF_Fortran_VARN)
    if (PnetCDF_Fortran_VARN)
      list (APPEND PnetCDF_Fortran_DEFINITIONS USE_PNETCDF_VARN)
      list (APPEND PnetCDF_Fortran_DEFINITIONS USE_PNETCDF_VARN_ON_READ)
    endif()  
endif ()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and 
# set PnetCDF_Fortran_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args (PnetCDF_Fortran DEFAULT_MSG
                                   PnetCDF_Fortran_LIBRARY PnetCDF_Fortran_INCLUDE_DIR)
mark_as_advanced (PnetCDF_Fortran_INCLUDE_DIR PnetCDF_Fortran_LIBRARY)

# HACK For bug in CMake v3.0:
set (PnetCDF_Fortran_FOUND ${PNETCDF_FORTRAN_FOUND})
