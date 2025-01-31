FIND_PATH(Parmetis_INCLUDE_DIR
          parmetis.h
          PATHS ${OMEGA_PARMETIS_ROOT}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

# If not defined, assume the METIS path isthe same as ParMETIS
if(NOT DEFINED OMEGA_METIS_ROOT)
   set(OMEGA_METIS_ROOT ${OMEGA_PARMETIS_ROOT})
endif()

FIND_PATH(Metis_INCLUDE_DIR
          metis.h
          PATHS ${OMEGA_METIS_ROOT}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

# Assume the GKlib path is the same as METIS if it is not defined.
# This library is not mandatory and is therefore optional.
if(NOT DEFINED OMEGA_GKLIB_ROOT)
   set(OMEGA_GKLIB_ROOT ${OMEGA_PARMETIS_ROOT})
endif()

FIND_PATH(GKlib_INCLUDE_DIR
          GKlib.h
	  PATHS ${OMEGA_GKLIB_ROOT}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

# Currently using static libraries, but retain the following for
# potential future use with shared libraries.
#IF (${OMEGA_PREFER_SHARED})
#  FIND_LIBRARY(Parmetis_LIBRARY
#               NAMES parmetis
#               HINTS ${Parmetis_INCLUDE_DIR}/../lib
#               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
#
#  FIND_LIBRARY(Metis_LIBRARY
#               NAMES metis
#               HINTS ${Metis_INCLUDE_DIR}/../lib
#               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
#
#  IF (${GKlib_INCLUDE_DIR})
#    FIND_LIBRARY(GKlib_LIBRARY
#                 NAMES GKlib
#	         HINTS ${GKlib_INCLUDE_DIR}/../lib
#	         NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
#  ENDIF ()
#
#ELSE ()

  FIND_LIBRARY(Parmetis_LIBRARY
               NAMES libparmetis.a
               HINTS ${Parmetis_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  FIND_LIBRARY(Metis_LIBRARY
               NAMES libmetis.a
               HINTS ${Metis_INCLUDE_DIR}/../lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  # In some installations, GKlib is optional
  IF (DEFINED GKlib_INCLUDE_DIR AND NOT GKlib_INCLUDE_DIR STREQUAL "")
    FIND_LIBRARY(GKlib_LIBRARY
                 NAMES libGKlib.a
                 HINTS ${GKlib_INCLUDE_DIR}/../lib
                 NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
  ENDIF ()

#ENDIF ()

SET(Parmetis_INCLUDE_DIRS)

IF (Parmetis_INCLUDE_DIR AND Parmetis_LIBRARY)

  SET(Parmetis_FOUND TRUE)
  LIST(APPEND Parmetis_INCLUDE_DIRS ${Parmetis_INCLUDE_DIR})

  MESSAGE(STATUS "Found Parmetis Library: ${Parmetis_LIBRARY}")
  MESSAGE(STATUS "Found Parmetis Include: ${Parmetis_INCLUDE_DIR}")
ELSE()
  SET(Parmetis_FOUND FALSE)
ENDIF()

IF (Metis_INCLUDE_DIR AND Metis_LIBRARY)

  LIST(APPEND Parmetis_INCLUDE_DIRS ${Metis_INCLUDE_DIR})
  SET(Metis_FOUND TRUE)

  MESSAGE(STATUS "Found Metis Library: ${Metis_LIBRARY}")
  MESSAGE(STATUS "Found Metis Include: ${Metis_INCLUDE_DIR}")
ELSE()
  SET(Metis_FOUND FALSE)
ENDIF()

IF (GKlib_INCLUDE_DIR AND GKlib_LIBRARY)

  LIST(APPEND Parmetis_INCLUDE_DIRS ${GKlib_INCLUDE_DIR})
  SET(GKlib_FOUND TRUE)

  MESSAGE(STATUS "Found GKlib Library: ${GKlib_LIBRARY}")
  MESSAGE(STATUS "Found GKlib Include: ${GKlib_INCLUDE_DIR}")
ELSE()
  SET(GKlib_FOUND FALSE)
ENDIF()

IF(NOT Parmetis_FOUND OR NOT Metis_FOUND)
  MESSAGE(STATUS "Did not find required library Parmetis.\n"
          "Please set location of Parmetis with -DOMEGA_PARMETIS_ROOT")
ENDIF()
