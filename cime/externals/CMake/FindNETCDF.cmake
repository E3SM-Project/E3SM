# - Try to find Netcdf
# Once done this will define
#  NETCDF_FOUND - System has Netcdf
#  NETCDF_INCLUDE_DIRS - The Netcdf include directories
# NETCDF_C_LIBRARIES - The C libraries needed to use Netcdf
# NETCDF_Fortran_LIBRARIES - The Fortran libraries needed to use Netcdf
#  NETCDF_LIBRARIES - All the libraries needed to use Netcdf
#  NETCDF_DEFINITIONS - Compiler switches required for using Netcdf

# If we weren't given a hint via a CMake variable, check the environment.
if(NOT NETCDF_DIR)
  set(NETCDF_DIR $ENV{NETCDF})
endif()

find_path(NETCDF_INCLUDE_DIR netcdf.h
          HINTS ${NETCDF_DIR}/include )

find_path(NETCDF_LIB_DIR NAMES libnetcdf.a libnetcdf.so
          HINTS ${NETCDF_DIR}/lib ${NETCDF_DIR}/lib64 )

find_path(NETCDF_FORTRAN_LIB_DIR NAMES libnetcdff.a libnetcdff.so
          HINTS ${NETCDF_DIR}/lib ${NETCDF_DIR}/lib64 )


find_file(NETCDF4_PAR_H netcdf_par.h
          HINTS ${NETCDF_INCLUDE_DIR}
          NO_DEFAULT_PATH )

#MESSAGE("PAR_H: ${NETCDF4_PAR_H}")
find_library(NETCDF_C_LIBRARY NAMES libnetcdf.a netcdf HINTS ${NETCDF_LIB_DIR})

if(NOT NETCDF_FORTRAN_LIB_DIR)
  MESSAGE(WARNING "Did not find netCDF Fortran library.")
else()
  find_library(NETCDF_Fortran_LIBRARY NAMES libnetcdff.a netcdff HINTS ${NETCDF_FORTRAN_LIB_DIR})
endif()
if(NOT NETCDF4_PAR_H)
  set(NETCDF4_PARALLEL "no")
  MESSAGE("NETCDF built without MPIIO")
else()
  set(NETCDF4_PARALLEL "yes")
  MESSAGE("NETCDF built with hdf5 MPIIO support")
endif()

set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR} )

FIND_PACKAGE(HDF5 COMPONENTS C HL CONFIG)

if(${HDF5_FOUND})
  MESSAGE(STATUS "Adding hdf5 libraries ")
 set(NETCDF_C_LIBRARY ${NETCDF_C_LIBRARY} ${HDF5_LIBRARIES}
                      ${SZIP_LIBRARIES} ${ZLIB_LIBRARIES})
endif()

# If netCDF was configured with DAP, it depends on libcurl.
find_program(NETCDF_NC_CONFIG nc-config HINTS ${NETCDF_INCLUDE_DIR}/../bin)
if(NETCDF_NC_CONFIG)
  execute_process(COMMAND ${NETCDF_NC_CONFIG} --has-dap
    OUTPUT_VARIABLE nc_config_output)
  if(nc_config_output MATCHES yes)
    find_package(CURL)
    if(CURL_FOUND)
      MESSAGE(STATUS "Adding curl libraries for netCDF DAP.")
      set(NETCDF_C_LIBRARY ${NETCDF_C_LIBRARY} ${CURL_LIBRARIES})
    else()
      MESSAGE(WARNING "netCDF DAP appears enabled, but libcurl was not found.")
    endif()
  endif()
endif()

set(NETCDF_LIBRARIES ${NETCDF_Fortran_LIBRARY} ${NETCDF_C_LIBRARY})

# Export variables so other projects can use them as well
#  ie. if pio is added with add_subdirectory
SET(NETCDF_INCLUDE_DIR ${NETCDF_INCLUDE_DIR} CACHE STRING "Location of NetCDF include files.")
SET(NETCDF_LIBRARIES ${NETCDF_LIBRARIES} CACHE STRING "Link line for NetCDF.")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE
# if all listed variables are TRUE
# (Note that the Fortran interface is not always a separate library, so
# don't require it to be found.)
find_package_handle_standard_args(NETCDF  DEFAULT_MSG NETCDF_LIBRARIES
                                  NETCDF_C_LIBRARY NETCDF_INCLUDE_DIR)

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_LIBRARIES NETCDF_C_LIBRARY NETCDF_Fortran_LIBRARY NETCDF4_PARALLEL )
