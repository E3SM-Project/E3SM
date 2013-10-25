


FIND_PATH(Scotch_INCLUDE_DIR
          scotch.h
          PATHS ${SCOTCH_DIR} ${Homme_SCOTCH_DIR}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

SET(Scotch_INCLUDE_DIRS ${Scotch_INCLUDE_DIR})

SET(Scotch_LIB_LIST
  ptscotch
  ptscotchparmetis
  ptscotcherr
  ptscotcherrexit
  scotch
  scotchmetis
  scotcherr
  scotcherrexit)
  
SET(Scotch_LIBRARIES)

FOREACH (Scotch_LIB ${Scotch_LIB_LIST})
  FIND_LIBRARY(Scotch_TMP_LIB
               NAMES ${Scotch_LIB}
               HINTS ${Scotch_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
  
  IF (Scotch_TMP_LIB)
    SET(Scotch_LIBRARIES ${Scotch_LIBRARIES} ${Scotch_TMP_LIB})
  ENDIF ()
  SET(Scotch_TMP_LIB)

ENDFOREACH ()

IF (Scotch_INCLUDE_DIR AND Scotch_LIBRARIES)
  SET(Scotch_FOUND TRUE)
  MESSAGE(STATUS "Found Scotch:")
  MESSAGE(STATUS "  Libraries: ${Scotch_LIBRARIES}")
  MESSAGE(STATUS "  Includes:  ${Scotch_INCLUDE_DIRS}")
ELSE()
  SET(Scotch_FOUND FALSE)
ENDIF()

IF(Scotch_FIND_REQUIRED AND NOT Scotch_FOUND)
  MESSAGE(STATUS "Did not find required library Scotch.\n"
          "Please set location of Scotch with -DSCOTCH_DIR")
ENDIF()
