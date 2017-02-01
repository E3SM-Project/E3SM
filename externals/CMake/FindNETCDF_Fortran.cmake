# - Try to find Netcdf
# Once done this will define
#  NETCDF_Fortran_FOUND - System has Netcdf
#  NETCDF_Fortran_INCLUDE_DIRS - The Netcdf include directories
#  NETCDF_Fortran_LIBRARY - The C libraries needed to use Netcdf
#  NETCDF_Fortran_LIBRARIES - All the libraries needed to use Netcdf
#  NETCDF_Fortran_DEFINITIONS - Compiler switches required for using Netcdf
#
include(LibFindMacros)

find_path(NETCDF_Fortran_INCLUDE_DIR
          NAMES netcdf.inc
          PATHS ${NETCDF_Fortran_PKGCONF_INCLUDE_DIRS}
          HINTS ${NETCDF_DIR}/include ${NETCDF_Fortran_DIR}/include)
libfind_package(NETCDF_Fortran NETCDF_C)
find_library(NETCDF_Fortran_LIBRARY
             NAMES libnetcdff.a netcdff
             PATHS ${NETCDF_Fortran_PKGCONF_LIBRARY_DIRS}
	     HINTS ${NETCDF_DIR}/lib ${NETCDF_Fortran_DIR}/lib)

set(NETCDF_Fortran_PROCESS_INCLUDES NETCDF_Fortran_INCLUDE_DIR)

set(NETCDF_Fortran_PROCESS_LIBS NETCDF_Fortran_LIBRARY)

libfind_process(NETCDF_Fortran)
