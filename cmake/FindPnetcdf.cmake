


find_path(Pnetcdf_INCLUDE_DIR 
          pnetcdf.h
          HINTS ${Homme_Hint_Paths}
          PATHS ${PNETCDF_DIR} ${Homme_PNETCDF_DIR}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

find_library(Pnetcdf_LIBRARY 
             NAMES libpnetcdf.a pnetcdf
             HINTS ${Pnetcdf_INCLUDE_DIR}/../lib
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

if (${Pnetcdf_INCLUDE_DIR} STREQUAL "Pnetcdf_INCLUDE_DIR-NOTFOUND")
  set(Pnetcdf_FOUND OFF)
else()
  set(Pnetcdf_LIBRARIES ${Pnetcdf_LIBRARY} )
  set(Pnetcdf_INCLUDE_DIRS ${Pnetcdf_INCLUDE_DIR} )
  set(Pnetcdf_FOUND ON)
  MESSAGE(STATUS "Found Pnetcdf:")
  MESSAGE(STATUS "  Libraries: ${Pnetcdf_LIBRARIES}")
  MESSAGE(STATUS "  Includes:  ${Pnetcdf_INCLUD_DIRS}")
endif()

IF(Pnetcdf_FIND_REQUIRED AND NOT Pnetcdf_FOUND)
  MESSAGE(FATAL_ERROR "Did not find required library Pnetcdf")
ENDIF()
