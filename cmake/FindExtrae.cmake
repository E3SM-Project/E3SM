
IF (DEFINED EXTRAE_LIB) 
  find_library(Extrae_LIBRARY 
               NAMES ${EXTRAE_LIB}
               HINTS ${EXTRAE_DIR} 
               PATH_SUFFIXES lib lib64
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

ELSE ()

  find_library(Extrae_LIBRARY 
               NAMES mpitracef
               HINTS ${EXTRAE_DIR} 
               PATH_SUFFIXES lib lib64
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

ENDIF ()

IF (Extrae_LIBRARY) 
  SET(Extrae_FOUND TRUE)
ELSE ()
  SET(Extrae_FOUND FALSE)
ENDIF ()


IF(Extrae_FIND_REQUIRED AND NOT Extrae_FOUND)
  MESSAGE(FATAL_ERROR "Did not find required library Extrae")
ELSE ()
  MESSAGE(STATUS "Found Extrae:")
  MESSAGE(STATUS "  Libraries: ${Extrae_LIBRARY}")
ENDIF()
