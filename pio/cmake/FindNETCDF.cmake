# - Try to find Netcdf
# Once done this will define
#  NETCDF_FOUND - System has Netcdf
#  NETCDF_INCLUDE_DIRS - The Netcdf include directories
#  NETCDF_LIBRARIES - The libraries needed to use Netcdf
#  NETCDF_DEFINITIONS - Compiler switches required for using Netcdf

find_path(NETCDF_INCLUDE_DIR netcdf.h
          HINTS ${NETCDF_DIR}/include )

find_library(NETCDF_FORTRAN_LIBRARY NAMES libnetcdff.a netcdff 
             HINTS ${NETCDF_DIR}/lib )

find_library(NETCDF_LIBRARY NAMES  libnetcdf.a netcdf
             HINTS ${NETCDF_DIR}/lib )

set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR} )

FIND_PACKAGE(HDF5 COMPONENTS C HL)

if(${HDF5_FOUND}) 
  MESSAGE(STATUS "Adding hdf5 libraries ")
 set(NETCDF_LIBRARIES ${NETCDF_FORTRAN_LIBRARY} ${NETCDF_LIBRARY} ${HDF5_LIBRARIES})  
else()
   set(NETCDF_LIBRARIES ${NETCDF_FORTRAN_LIBRARY} ${NETCDF_LIBRARY}  )  
endif()

TRY_COMPILE(NETCDF4_PARALLEL ${CMAKE_CURRENT_BINARY_DIR}/tryNC4
                          ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TryNC4.f90
			  COMPILE_DEFINITIONS -I${NETCDF_INCLUDE_DIR}
			  CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${NETCDF_LIBRARIES}"
                           OUTPUT_VARIABLE NC4_OUT)
if(NETCDF4_PARALLEL)
  MESSAGE(STATUS "netcdf4 built with parallel support ")
else()
  MESSAGE("netcdf4 not built with parallel support ${NC4_OUT}")
endif()


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NETCDF  DEFAULT_MSG 
                                  NETCDF_LIBRARY NETCDF_INCLUDE_DIR)

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_LIBRARY NETCDF4_PARALLEL)