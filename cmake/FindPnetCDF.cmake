# - Try to find PnetCDF
#
# Once done, this will define (for each supported language)
#
#  PnetCDF_<lang>_FOUND        - system has PnetCDF
#  PnetCDF_<lang>_INCLUDE_DIRS - the PnetCDF include directories
#  PnetCDF_<lang>_LIBRARIES    - link these to use PnetCDF

include(LibFindMacros)

#=== C Language Package ===

# Dependencies
libfind_package(PnetCDF_C MPI)

# Include dir
set (pnetcdf_inc_search_paths)
if (PnetCDF_C_INCLUDE_DIR)
    list (APPEND pnetcdf_inc_search_paths ${PnetCDF_C_INCLUDE_DIR})
endif ()
if (PnetCDF_DIR)
    list (APPEND pnetcdf_inc_search_paths ${PnetCDF_DIR}/include)
endif ()
if (ENV{PNETCDF})
    list (APPEND pnetcdf_inc_search_paths $ENV{PNETCDF}/include)
endif ()
find_path(PnetCDF_C_INCLUDE_DIR
  NAMES pnetcdf.h
  PATHS ${pnetcdf_inc_search_paths})
unset (pnetcdf_inc_search_paths)

# Library itself
set (pnetcdf_lib_search_paths)
if (PnetCDF_C_LIBRARY)
    get_filename_component (pnetcdf_lib_dir ${PnetCDF_C_LIBRARY} DIRECTORY)
    list (APPEND pnetcdf_lib_search_paths ${pnetcdf_lib_dir})
endif ()
if (PnetCDF_DIR)
    list (APPEND pnetcdf_lib_search_paths ${PnetCDF_DIR}/lib)
endif ()
if (ENV{PNETCDF})
    list (APPEND pnetcdf_lib_search_paths $ENV{PNETCDF}/lib)
endif ()
find_library(PnetCDF_C_LIBRARY
  NAMES pnetcdf libpnetcdf.a
  PATHS ${pnetcdf_lib_search_paths})
unset (pnetcdf_lib_search_paths)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(PnetCDF_C_PROCESS_INCLUDES PnetCDF_C_INCLUDE_DIR PnetCDF_C_INCLUDE_DIRS)
set(PnetCDF_C_PROCESS_LIBS PnetCDF_C_LIBRARY PnetCDF_C_LIBRARIES)
libfind_process(PnetCDF_C)

#=== C++ Language Package ===

# Dependencies
libfind_package(PnetCDF_CXX MPI)

# Include dir
set (pnetcdf_inc_search_paths)
if (PnetCDF_CXX_INCLUDE_DIR)
    list (APPEND pnetcdf_inc_search_paths ${PnetCDF_CXX_INCLUDE_DIR})
endif ()
if (PnetCDF_DIR)
    list (APPEND pnetcdf_inc_search_paths ${PnetCDF_DIR}/include)
endif ()
if (ENV{PNETCDF})
    list (APPEND pnetcdf_inc_search_paths $ENV{PNETCDF}/include)
endif ()
find_path(PnetCDF_CXX_INCLUDE_DIR
  NAMES pnetcdf
  PATHS ${pnetcdf_inc_search_paths})
unset (pnetcdf_inc_search_paths)

# Library itself
set (pnetcdf_lib_search_paths)
if (PnetCDF_CXX_LIBRARY)
    get_filename_component (pnetcdf_lib_dir ${PnetCDF_CXX_LIBRARY} DIRECTORY)
    list (APPEND pnetcdf_lib_search_paths ${pnetcdf_lib_dir})
endif ()
if (PnetCDF_DIR)
    list (APPEND pnetcdf_lib_search_paths ${PnetCDF_DIR}/lib)
endif ()
if (ENV{PNETCDF})
    list (APPEND pnetcdf_lib_search_paths $ENV{PNETCDF}/lib)
endif ()
find_library(PnetCDF_CXX_LIBRARY
  NAMES pnetcdf libpnetcdf.a
  PATHS ${pnetcdf_lib_search_paths})
unset (pnetcdf_lib_search_paths)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(PnetCDF_CXX_PROCESS_INCLUDES PnetCDF_CXX_INCLUDE_DIR PnetCDF_CXX_INCLUDE_DIRS)
set(PnetCDF_CXX_PROCESS_LIBS PnetCDF_CXX_LIBRARY PnetCDF_CXX_LIBRARIES)
libfind_process(PnetCDF_CXX)

#=== Fortran Language Package ===

# Dependencies
libfind_package(PnetCDF_Fortran MPI)

# Include dir
set (pnetcdf_inc_search_paths)
if (PnetCDF_Fortran_INCLUDE_DIR)
    list (APPEND pnetcdf_inc_search_paths ${PnetCDF_Fortran_INCLUDE_DIR})
endif ()
if (PnetCDF_DIR)
    list (APPEND pnetcdf_inc_search_paths ${PnetCDF_DIR}/include)
endif ()
if (ENV{PNETCDF})
    list (APPEND pnetcdf_inc_search_paths $ENV{PNETCDF}/include)
endif ()
find_path(PnetCDF_Fortran_INCLUDE_DIR
  NAMES pnetcdf.mod pnetcdf.inc
  PATHS ${pnetcdf_inc_search_paths})
unset (pnetcdf_inc_search_paths)

# Library itself
set (pnetcdf_lib_search_paths)
if (PnetCDF_Fortran_LIBRARY)
    get_filename_component (pnetcdf_lib_dir ${PnetCDF_Fortran_LIBRARY} DIRECTORY)
    list (APPEND pnetcdf_lib_search_paths ${pnetcdf_lib_dir})
endif ()
if (PnetCDF_DIR)
    list (APPEND pnetcdf_lib_search_paths ${PnetCDF_DIR}/lib)
endif ()
if (ENV{PNETCDF})
    list (APPEND pnetcdf_lib_search_paths $ENV{PNETCDF}/lib)
endif ()
find_library(PnetCDF_Fortran_LIBRARY
  NAMES pnetcdf libpnetcdf.a
  PATHS ${pnetcdf_lib_search_paths})
unset (pnetcdf_lib_search_paths)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(PnetCDF_Fortran_PROCESS_INCLUDES PnetCDF_Fortran_INCLUDE_DIR PnetCDF_Fortran_INCLUDE_DIRS)
set(PnetCDF_Fortran_PROCESS_LIBS PnetCDF_Fortran_LIBRARY PnetCDF_Fortran_LIBRARIES)
libfind_process(PnetCDF_Fortran)

