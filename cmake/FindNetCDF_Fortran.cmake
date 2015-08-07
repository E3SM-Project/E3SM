# - Try to find NetCDF-Fortran
#
# This can be controlled by setting the NetCDF_DIR or NetCDF_Fortran_DIR
# CMake variables, or the NETCDF environment variable.  
#
# Once done, this will define:
#
#   NetCDF_Fortran_IS_SHARED    (BOOL) - Whether library is shared/dynamic
#   NetCDF_Fortran_FOUND        (BOOL) - system has NetCDF
#   NetCDF_Fortran_INCLUDE_DIR  (PATH) - Location of the Fortran include files
#   NetCDF_Fortran_INCLUDE_DIRS (LIST) - the NetCDF include directories
#   NetCDF_Fortran_LIBRARY      (FILE) - Path to the Fortran library file
#   NetCDF_Fortran_LIBRARIES    (LIST) - link these to use NetCDF
#   NetCDF_Fortran_DEFINITIONS  (LIST) - preprocessor macros to use with NetCDF
#   NetCDF_Fortran_OPTIONS      (LIST) - compiler options to use NetCDF

# Determine include dir search order
set (NetCDF_Fortran_INCLUDE_HINTS)
if (NetCDF_Fortran_INCLUDE_DIR)
    list (APPEND NetCDF_Fortran_INCLUDE_HINTS ${NetCDF_Fortran_INCLUDE_DIR})
endif ()
if (NetCDF_Fortran_DIR)
    list (APPEND NetCDF_Fortran_INCLUDE_HINTS ${NetCDF_Fortran_DIR}/include)
endif ()
if (NetCDF_DIR)
    list (APPEND NetCDF_Fortran_INCLUDE_HINTS ${NetCDF_DIR}/include)
endif ()
if (DEFINED ENV{NETCDF})
    list (APPEND NetCDF_Fortran_INCLUDE_HINTS $ENV{NETCDF}/include)
endif ()

# Search for include file
find_path (NetCDF_Fortran_INCLUDE_DIR
           NAMES netcdf.inc netcdf.mod
           HINTS ${NetCDF_Fortran_INCLUDE_HINTS})
           
# Unset include search variables
unset (NetCDF_Fortran_INCLUDE_HINTS)

# Determine library dir search order
set (NetCDF_Fortran_LIBRARY_HINTS)
if (NetCDF_Fortran_LIBRARY)
    get_filename_component (NetCDF_Fortran_library_path ${NetCDF_Fortran_LIBRARY} PATH)
    list (APPEND NetCDF_Fortran_LIBRARY_HINTS ${NetCDF_Fortran_library_path})
    unset (NetCDF_Fortran_library_path)
endif ()
if (NetCDF_Fortran_DIR)
    list (APPEND NetCDF_Fortran_LIBRARY_HINTS ${NetCDF_Fortran_DIR}/lib)
endif ()
if (NetCDF_DIR)
    list (APPEND NetCDF_Fortran_LIBRARY_HINTS ${NetCDF_DIR}/lib)
endif ()
if (DEFINED ENV{NETCDF})
    list (APPEND NetCDF_Fortran_LIBRARY_HINTS $ENV{NETCDF}/lib)
endif ()

# Search for library file
set (NetCDF_Fortran_IS_SHARED FALSE)
include (LibFindLibraryMacros)
find_shared_library (NetCDF_Fortran_LIBRARY
                     NAMES netcdff
                     HINTS ${NetCDF_Fortran_LIBRARY_HINTS}
                     NO_DEFAULT_PATH)
if (NetCDF_Fortran_LIBRARY)
    set (NetCDF_Fortran_IS_SHARED TRUE)
else ()
    find_static_library (NetCDF_Fortran_LIBRARY
                         NAMES netcdff
                         HINTS ${NetCDF_Fortran_LIBRARY_HINTS}
                         NO_DEFAULT_PATH)
endif ()

# Unset include search variables
unset (NetCDF_Fortran_LIBRARY_HINTS)

# Set return variables
set (NetCDF_Fortran_INCLUDE_DIRS ${NetCDF_Fortran_INCLUDE_DIR})
set (NetCDF_Fortran_LIBRARIES ${NetCDF_Fortran_LIBRARY})
set (NetCDF_Fortran_DEFINITIONS)
set (NetCDF_Fortran_OPTIONS)

# If static, look for dependencies
if (NOT NetCDF_Fortran_IS_SHARED)

    # DEPENDENCY: HDF5
    find_package (HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)
    if (HDF5_FOUND)
        list (APPEND NetCDF_Fortran_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
        list (APPEND NetCDF_Fortran_LIBRARIES ${HDF5_Fortran_LIBRARIES}
                                              ${HDF5_Fortran_HL_LIBRARIES})
    endif ()

    # DEPENDENCY: CURL
    find_package (CURL REQUIRED)
    if (CURL_FOUND)
        list (APPEND NetCDF_Fortran_INCLUDE_DIRS ${CURL_INCLUDE_DIRS})
        list (APPEND NetCDF_Fortran_LIBRARIES ${CURL_LIBRARIES})
    endif ()
     
endif ()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and 
# set NetCDF_Fortran_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args (NetCDF_Fortran DEFAULT_MSG
                                   NetCDF_Fortran_LIBRARY NetCDF_Fortran_INCLUDE_DIR)
mark_as_advanced (NetCDF_Fortran_INCLUDE_DIR NetCDF_Fortran_LIBRARY)

# HACK For bug in CMake v3.0:
set (NetCDF_Fortran_FOUND ${NETCDF_FORTRAN_FOUND})
