PROJECT(PIO C Fortran)
ENABLE_LANGUAGE(Fortran)
INCLUDE(FortranCInterface)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
INCLUDE(SetHostPaths)
# Netcdf is required
FIND_PACKAGE(Netcdf REQUIRED)
# PNetcdf is optional but used by default
OPTION(WITH_PNETCDF "Whether to build with PNetcdf" TRUE)
IF (${WITH_PNETCDF})
  FIND_PACKAGE(Pnetcdf REQUIRED)
ELSE ()
  MESSAGE(WARNING "Warning: Not building with PNetcdf - cannot run all regression tests.")
ENDIF ()

SET (GENF90_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

ADD_SUBDIRECTORY(pio)
ADD_SUBDIRECTORY(testpio)
