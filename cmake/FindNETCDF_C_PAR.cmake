# - Try to find NETCDF_C_PAR
#
# Once done, this will define
#
#  NETCDF_C_PAR_FOUND        - system has NETCDF_C_PAR
#  NETCDF_C_PAR_INCLUDE_DIRS - the NETCDF_C_PAR include directories
#  NETCDF_C_PAR_LIBRARIES    - link these to use NETCDF_C_PAR

include(LibFindMacros)

# Dependencies
libfind_package(NETCDF_C_PAR NETCDF_C)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(NETCDF_C_PAR_PKGCONF NETCDF_C_PAR)

# Include dir
find_path(NETCDF_C_PAR_INCLUDE_DIR
  NAMES netcdf_par.h
  PATHS ${NETCDF_C_PAR_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(NETCDF_C_PAR_LIBRARY
  NAMES libnetcdf.a netcdf
  PATHS ${NETCDF_C_PAR_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(NETCDF_C_PAR_PROCESS_INCLUDES NETCDF_C_PAR_INCLUDE_DIR NETCDF_C_PAR_INCLUDE_DIRS)
set(NETCDF_C_PAR_PROCESS_LIBS NETCDF_C_PAR_LIBRARY NETCDF_C_PAR_LIBRARIES)
libfind_process(NETCDF_C_PAR)