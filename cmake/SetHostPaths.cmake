
# Determine the host name
EXECUTE_PROCESS(COMMAND uname -n
  RESULT_VARIABLE Homme_result
  OUTPUT_VARIABLE Homme_output
  ERROR_VARIABLE Homme_error
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

IF (Homme_result EQUAL 0 AND Homme_error STREQUAL "")
  SET(Homme_Raw_Hostname ${Homme_output})
  #MESSAGE(STATUS "Raw Hostname = ${Homme_Raw_Hostname}")
ELSE ()
  MESSAGE(STATUS "Hostname could not be determined")
ENDIF()


EXECUTE_PROCESS(COMMAND uname -s
  RESULT_VARIABLE Homme_result
  OUTPUT_VARIABLE Homme_output
  ERROR_VARIABLE Homme_error
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

IF (Homme_result EQUAL 0 AND Homme_error STREQUAL "")
  SET(Homme_OS ${Homme_output})
ELSE ()
  MESSAGE(STATUS "OS type could not be determined")
ENDIF()

EXECUTE_PROCESS(COMMAND whoami
  RESULT_VARIABLE Homme_result
  OUTPUT_VARIABLE Homme_output
  ERROR_VARIABLE Homme_error
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

IF (Homme_result EQUAL 0 AND Homme_error STREQUAL "")
  SET(Homme_Username ${Homme_output})
  #MESSAGE(STATUS "Homme_Username = ${Homme_Username}")
ELSE ()
  MESSAGE(STATUS "Username type could not be determined")
ENDIF()

# Now parse the host name to see if it matches yslogin*
# MATCHES Does REGEX where "." is a wildcard
IF (${Homme_Raw_Hostname} MATCHES "yslogin.")
  SET(Homme_Hostname "Yellowstone")
  SET(Homme_Registered_Host TRUE)
ELSEIF (${Homme_Raw_Hostname} MATCHES "titan-ext." OR ${Homme_Raw_Hostname} MATCHES "titan-login.")
  SET(Homme_Hostname "Titan")
  SET(Homme_Registered_Host TRUE)
ELSE ()
  SET(Homme_Hostname ${Homme_Raw_Hostname})
  SET(Homme_Registered_Host FALSE)
ENDIF()

# Back up cached variables 
get_cmake_property(CACHE_VARS CACHE_VARIABLES)
foreach(CACHE_VAR ${CACHE_VARS})
  get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
  if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line.")
    set(${CACHE_VAR}_BACK ${${CACHE_VAR}})
    set(ORIG_LIST ${ORIG_LIST};${CACHE_VAR})
    set(BACKUP_LIST ${BACKUP_LIST};"${CACHE_VAR}_BACK")
  endif()
endforeach()

IF (Homme_Registered_Host)
  # Try to read the system variables from the machinesFiles
  MESSAGE(STATUS "Registered Host, reading machineFile")
  STRING(TOLOWER ${Homme_Hostname} MACHINEFILE)
  INCLUDE(machineFiles/${MACHINEFILE} OPTIONAL RESULT_VARIABLE MACH_FILE_FOUND)
  IF (NOT MACH_FILE_FOUND)

    SET(MACHINEFILE_COMP ${MACHINEFILE}${CMAKE_Fortran_COMPILER_ID})
    
    INCLUDE(machineFiles/${MACHINEFILE_COMP} OPTIONAL RESULT_VARIABLE MACH_FILE_FOUND)

    IF (NOT MACH_FILE_FOUND)
      MESSAGE(STATUS "Registered Host ${Homme_Hostname}, could not find a machine file")
      MESSAGE(STATUS "Looked at ${MACHINEFILE} and ${MACHINEFILE_COMP}")
      
    ELSE ()
      MESSAGE(STATUS "Reading machine specific info from ${MACHINEFILE_COMP}")
    ENDIF ()
  ELSE ()
    MESSAGE(STATUS "Reading machine specific info from ${MACHINEFILE}")
  ENDIF ()
ELSE ()
  # Set some possible paths for the OS type
  MESSAGE(STATUS "Host not registered setting hints ")
ENDIF ()

# Now set back the command line specified variables
foreach(CACHE_VAR ${ORIG_LIST})
  SET(${CACHE_VAR} ${${CACHE_VAR}_BACK})
endforeach()
