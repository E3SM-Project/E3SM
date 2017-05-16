# - Try to find Netcdf
# Once done this will define
#  NETCDF_C_FOUND - System has Netcdf
#  NETCDF_C_INCLUDE_DIRS - The Netcdf include directories
#  NETCDF_C_LIBRARY - The C libraries needed to use Netcdf
#  NETCDF_C_LIBRARIES - All the libraries needed to use Netcdf
#  NETCDF_C_DEFINITIONS - Compiler switches required for using Netcdf
#
include(LibFindMacros)
set(NETCDF_C_PARALLEL FALSE)
find_path(NETCDF_C_INCLUDE_DIR
          NAMES netcdf.h
          PATHS ${NETCDF_C_PKGCONF_INCLUDE_DIRS}
          HINTS ${NETCDF_DIR}/include ${NETCDF_C_DIR}/include)

# See if netcdf includes parallel support
find_path(NETCDF_C_PAR_INCLUDE_DIR
          NAMES netcdf_par.h
          PATHS ${NETCDF_C_PKGCONF_INCLUDE_DIRS}
          HINTS ${NETCDF_DIR}/include ${NETCDF_C_DIR}/include)

set(NETCDF_C_DEFINITIONS "-D_NETCDF")
if(${NETCDF_C_PAR_INCLUDE_DIR} STREQUAL "NETCDF_C_PAR_INCLUDE_DIR-NOTFOUND")
  MESSAGE("Netcdf library does not appear to have parallel IO support")
else()
  MESSAGE("Netcdf library includes HDF5 parallel support")
  LIST(APPEND NETCDF_C_DEFINITIONS "-D_NETCDF4")
endif()
find_library(NETCDF_C_LIBRARY
             NAMES libnetcdf.a netcdf
             PATHS ${NETCDF_C_PKGCONF_LIBRARY_DIRS}
	     HINTS ${NETCDF_DIR}/lib ${NETCDF_C_DIR}/lib)

set(NETCDF_C_PROCESS_INCLUDES NETCDF_C_INCLUDE_DIR)

set(NETCDF_C_PROCESS_LIBS NETCDF_C_LIBRARY)

libfind_process(NETCDF_C)
