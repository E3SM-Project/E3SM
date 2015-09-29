


FIND_PATH(Parmetis_INCLUDE_DIR 
          parmetis.h
          PATHS ${PARMETIS_DIR} ${Homme_PARMETIS_DIR}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

IF (${PREFER_SHARED})
  FIND_LIBRARY(Parmetis_LIBRARY
               NAMES parmetis
               HINTS ${Parmetis_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  FIND_LIBRARY(Metis_LIBRARY
               NAMES metis
               HINTS ${Parmetis_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

ELSE ()
  FIND_LIBRARY(Parmetis_LIBRARY 
               NAMES libparmetis.a 
               HINTS ${Parmetis_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  FIND_LIBRARY(Metis_LIBRARY 
               NAMES libmetis.a 
               HINTS ${Parmetis_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

ENDIF ()

IF (Parmetis_INCLUDE_DIR AND Parmetis_LIBRARY AND Metis_LIBRARY)
  SET(Parmetis_LIBRARIES ${Parmetis_LIBRARY} ${Metis_LIBRARY})
  SET(Parmetis_INCLUDE_DIRS ${Parmetis_INCLUDE_DIR} )
  SET(Parmetis_FOUND TRUE)
  MESSAGE(STATUS "Found Parmetis:")
  MESSAGE(STATUS "  Libraries: ${Parmetis_LIBRARIES}")
  MESSAGE(STATUS "  Includes:  ${Parmetis_INCLUDE_DIRS}")
ELSE()
  SET(Parmetis_FOUND FALSE)
ENDIF()

IF(Parmetis_FIND_REQUIRED AND NOT Parmetis_FOUND)
  MESSAGE(STATUS "Did not find required library Parmetis.\n"
          "Please set location of Parmetis with -DPARMETIS_DIR")
ENDIF()
