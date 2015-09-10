


FIND_PATH(Pnetcdf_INCLUDE_DIR 
          pnetcdf.h
          HINTS ${Homme_Hint_Paths}
          PATHS ${PNETCDF_DIR} ${Homme_PNETCDF_DIR}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

IF (${PREFER_SHARED})
  FIND_LIBRARY(Pnetcdf_LIBRARY 
               NAMES pnetcdf
               HINTS ${Pnetcdf_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
ELSE ()
  FIND_LIBRARY(Pnetcdf_LIBRARY 
               NAMES libpnetcdf.a pnetcdf
               HINTS ${Pnetcdf_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
ENDIF ()

IF (Pnetcdf_INCLUDE_DIR AND Pnetcdf_LIBRARY)
  SET(Pnetcdf_LIBRARIES ${Pnetcdf_LIBRARY} )
  SET(Pnetcdf_INCLUDE_DIRS ${Pnetcdf_INCLUDE_DIR} )
  SET(Pnetcdf_FOUND ON)
  MESSAGE(STATUS "Found Pnetcdf:")
  MESSAGE(STATUS "  Libraries: ${Pnetcdf_LIBRARIES}")
  MESSAGE(STATUS "  Includes:  ${Pnetcdf_INCLUDE_DIRS}")
ELSE()
  SET(Pnetcdf_FOUND OFF)
ENDIF()

IF(Pnetcdf_FIND_REQUIRED AND NOT Pnetcdf_FOUND)
  MESSAGE(FATAL_ERROR "Did not find required library Pnetcdf.\n"
          "Please set location of Pnetcdf with -DPNETCDF_DIR")
ENDIF()
