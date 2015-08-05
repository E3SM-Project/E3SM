# - Try to find NetCDF
#
# Once done, this will define (for each supported language)
#
#  NetCDF_<lang>_FOUND        - system has NetCDF
#  NetCDF_<lang>_INCLUDE_DIRS - the NetCDF include directories
#  NetCDF_<lang>_LIBRARIES    - link these to use NetCDF

include(LibFindMacros)

#=== C Language Package ===

# Dependencies
libfind_package(NetCDF_C MPI)

# Include dir
set (netcdf_inc_search_paths)
if (NetCDF_C_INCLUDE_DIR)
    list (APPEND netcdf_inc_search_paths ${NetCDF_C_INCLUDE_DIR})
endif ()
if (NetCDF_DIR)
    list (APPEND netcdf_inc_search_paths ${NetCDF_DIR}/include)
endif ()
if (ENV{NETCDF})
    list (APPEND netcdf_inc_search_paths $ENV{NETCDF}/include)
endif ()
find_path(NetCDF_C_INCLUDE_DIR
    NAMES netcdf_par.h netcdf.h
    PATHS ${netcdf_inc_search_paths})
unset (netcdf_inc_search_paths)

# Library itself
set (netcdf_lib_search_paths)
if (NetCDF_C_LIBRARY)
    get_filename_component (netcdf_lib_dir ${NetCDF_C_LIBRARY} DIRECTORY)
    list (APPEND netcdf_lib_search_paths ${netcdf_lib_dir})
endif ()
if (NetCDF_DIR)
    list (APPEND netcdf_lib_search_paths ${NetCDF_DIR}/lib)
endif ()
if (ENV{NETCDF})
    list (APPEND netcdf_lib_search_paths $ENV{NETCDF}/lib)
endif ()
find_library(NetCDF_C_LIBRARY
    NAMES netcdf libnetcdf.a
    PATHS ${netcdf_lib_search_paths})
unset (netcdf_lib_search_paths)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(NetCDF_C_PROCESS_INCLUDES NetCDF_C_INCLUDE_DIR NetCDF_C_INCLUDE_DIRS)
set(NetCDF_C_PROCESS_LIBS NetCDF_C_LIBRARY NetCDF_C_LIBRARIES)
libfind_process(NetCDF_C)

#=== C++ Language Package ===

# Dependencies
libfind_package(NetCDF_CXX MPI)

# Include dir
set (netcdf_inc_search_paths)
if (NetCDF_CXX_INCLUDE_DIR)
    list (APPEND netcdf_inc_search_paths ${NetCDF_CXX_INCLUDE_DIR})
endif ()
if (NetCDF_DIR)
    list (APPEND netcdf_inc_search_paths ${NetCDF_DIR}/include)
endif ()
if (ENV{NETCDF})
    list (APPEND netcdf_inc_search_paths $ENV{NETCDF}/include)
endif ()
find_path(NetCDF_CXX_INCLUDE_DIR
    NAMES netcdf
    PATHS ${netcdf_inc_search_paths})
unset (netcdf_inc_search_paths)

# Library itself
set (netcdf_lib_search_paths)
if (NetCDF_CXX_LIBRARY)
    get_filename_component (netcdf_lib_dir ${NetCDF_CXX_LIBRARY} DIRECTORY)
    list (APPEND netcdf_lib_search_paths ${netcdf_lib_dir})
endif ()
if (NetCDF_DIR)
    list (APPEND netcdf_lib_search_paths ${NetCDF_DIR}/lib)
endif ()
if (ENV{NETCDF})
    list (APPEND netcdf_lib_search_paths $ENV{NETCDF}/lib)
endif ()
find_library(NetCDF_CXX_LIBRARY
    NAMES netcdf_c++ libnetcdf_c++.a
    PATHS ${netcdf_lib_search_paths})
unset (netcdf_lib_search_paths)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(NetCDF_CXX_PROCESS_INCLUDES NetCDF_CXX_INCLUDE_DIR NetCDF_CXX_INCLUDE_DIRS)
set(NetCDF_CXX_PROCESS_LIBS NetCDF_CXX_LIBRARY NetCDF_CXX_LIBRARIES)
libfind_process(NetCDF_CXX)

#=== Fortran Language Package ===

# Dependencies
libfind_package(NetCDF_Fortran MPI)

# Include dir
set (netcdf_inc_search_paths)
if (NetCDF_Fortran_INCLUDE_DIR)
    list (APPEND netcdf_inc_search_paths ${NetCDF_Fortran_INCLUDE_DIR})
endif ()
if (NetCDF_DIR)
    list (APPEND netcdf_inc_search_paths ${NetCDF_DIR}/include)
endif ()
if (ENV{NETCDF})
    list (APPEND netcdf_inc_search_paths $ENV{NETCDF}/include)
endif ()
find_path(NetCDF_Fortran_INCLUDE_DIR
    NAMES netcdf.mod netcdf.inc
    PATHS ${netcdf_inc_search_paths})
unset (netcdf_inc_search_paths)

# Library itself
set (netcdf_lib_search_paths)
if (NetCDF_Fortran_LIBRARY)
    get_filename_component (netcdf_lib_dir ${NetCDF_Fortran_LIBRARY} DIRECTORY)
    list (APPEND netcdf_lib_search_paths ${netcdf_lib_dir})
endif ()
if (NetCDF_DIR)
    list (APPEND netcdf_lib_search_paths ${NetCDF_DIR}/lib)
endif ()
if (ENV{NETCDF})
    list (APPEND netcdf_lib_search_paths $ENV{NETCDF}/lib)
endif ()
find_library(NetCDF_Fortran_LIBRARY
    NAMES netcdff libnetcdff.a
    PATHS ${netcdf_lib_search_paths})
unset (netcdf_lib_search_paths)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(NetCDF_Fortran_PROCESS_INCLUDES NetCDF_Fortran_INCLUDE_DIR NetCDF_Fortran_INCLUDE_DIRS)
set(NetCDF_Fortran_PROCESS_LIBS NetCDF_Fortran_LIBRARY NetCDF_Fortran_LIBRARIES)
libfind_process(NetCDF_Fortran)

