


FIND_PATH(Zoltan_INCLUDE_DIR 
          zoltan.h
          PATHS ${ZOLTAN_DIR} ${Homme_ZOLTAN_DIR}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

IF (${PREFER_SHARED})
  FIND_LIBRARY(Zoltan_LIBRARY
               NAMES zoltan
               HINTS ${Zoltan_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
ELSE ()
  FIND_LIBRARY(Zoltan_LIBRARY 
               NAMES libzoltan.a zoltan
               HINTS ${Zoltan_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
ENDIF ()

IF (Zoltan_INCLUDE_DIR AND Zoltan_LIBRARY)
  SET(Zoltan_LIBRARIES ${Zoltan_LIBRARY} )
  SET(Zoltan_INCLUDE_DIRS ${Zoltan_INCLUDE_DIR} )
  SET(Zoltan_FOUND TRUE)
  MESSAGE(STATUS "Found Zoltan:")
  MESSAGE(STATUS "  Libraries: ${Zoltan_LIBRARIES}")
  MESSAGE(STATUS "  Includes:  ${Zoltan_INCLUDE_DIRS}")
ELSE()
  SET(Zoltan_FOUND FALSE)
ENDIF()

IF(Zoltan_FIND_REQUIRED AND NOT Zoltan_FOUND)
  MESSAGE(STATUS "Did not find required library Zoltan.\n"
          "Please set location of Zoltan with -DZOLTAN_DIR")
ENDIF()
