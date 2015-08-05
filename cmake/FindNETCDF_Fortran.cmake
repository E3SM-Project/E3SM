# - Try to find NETCDF_Fortran
#
# Once done, this will define
#
#  NETCDF_Fortran_FOUND        - system has NETCDF_Fortran
#  NETCDF_Fortran_INCLUDE_DIRS - the NETCDF_Fortran include directories
#  NETCDF_Fortran_LIBRARIES    - link these to use NETCDF_Fortran

include(LibFindMacros)

# Dependencies
libfind_package(NETCDF_Fortran NETCDF_C)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(NETCDF_Fortran_PKGCONF NETCDF_Fortran)

# Include dir
find_path(NETCDF_Fortran_INCLUDE_DIR
  NAMES netcdf.inc
  PATHS ${NETCDF_Fortran_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(NETCDF_Fortran_LIBRARY
  NAMES libnetcdff.a netcdff
  PATHS ${NETCDF_Fortran_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(NETCDF_Fortran_PROCESS_INCLUDES NETCDF_Fortran_INCLUDE_DIR NETCDF_Fortran_INCLUDE_DIRS)
set(NETCDF_Fortran_PROCESS_LIBS NETCDF_Fortran_LIBRARY NETCDF_Fortran_LIBRARIES)
libfind_process(NETCDF_Fortran)