# - Try to find NETCDF_C
#
# Once done, this will define
#
#  NETCDF_C_FOUND        - system has NETCDF_C
#  NETCDF_C_INCLUDE_DIRS - the NETCDF_C include directories
#  NETCDF_C_LIBRARIES    - link these to use NETCDF_C

include(LibFindMacros)

# Dependencies
libfind_package(NETCDF_C HDF5)
libfind_package(NETCDF_C ZLIB)
libfind_package(NETCDF_C CURL)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(NETCDF_C_PKGCONF NETCDF_C)

# Include dir
find_path(NETCDF_C_INCLUDE_DIR
  NAMES netcdf.h
  PATHS ${NETCDF_C_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(NETCDF_C_LIBRARY
  NAMES libnetcdf.a netcdf
  PATHS ${NETCDF_C_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(NETCDF_C_PROCESS_INCLUDES NETCDF_C_INCLUDE_DIR NETCDF_C_INCLUDE_DIRS)
set(NETCDF_C_PROCESS_LIBS NETCDF_C_LIBRARY NETCDF_C_LIBRARIES)
libfind_process(NETCDF_C)