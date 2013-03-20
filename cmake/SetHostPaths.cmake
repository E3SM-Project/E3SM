# This file is used primarily to provide macro definitions and set the default locations for dependencies on supported systems
INCLUDE(DeterminePaths)

# Determine the host name
EXECUTE_PROCESS(COMMAND uname -n
  RESULT_VARIABLE Homme_result
  OUTPUT_VARIABLE Homme_output
  ERROR_VARIABLE Homme_error
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

IF (Homme_result EQUAL 0 AND Homme_error STREQUAL "")
  SET(Homme_Raw_Hostname ${Homme_output})
  MESSAGE(STATUS "Raw Hostname = ${Homme_Raw_Hostname}")
ELSE ()
  MESSAGE(FATAL_ERROR "Hostname could not be determined")
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
  MESSAGE(FATAL_ERROR "OS type could not be determined")
ENDIF()

EXECUTE_PROCESS(COMMAND whoami
  RESULT_VARIABLE Homme_result
  OUTPUT_VARIABLE Homme_output
  ERROR_VARIABLE Homme_error
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

IF (Homme_result EQUAL 0 AND Homme_error STREQUAL "")
  SET(Homme_Username ${Homme_output})
  MESSAGE(STATUS "Homme_Username = ${Homme_Username}")
ELSE ()
  MESSAGE(FATAL_ERROR "Username type could not be determined")
ENDIF()

# Now parse the host name to see if it matches yslogin*
# MATCHES Does REGEX where "." is a wildcard
IF (${Homme_Raw_Hostname} MATCHES "yslogin.")
  SET(Homme_Hostname "Yellowstone")
  SET(Homme_Registered_Host TRUE)
ELSEIF (${Homme_Raw_Hostname} MATCHES "titan-ext.")
  SET(Homme_Hostname "Titan")
  SET(Homme_Registered_Host FALSE)
ELSEIF (${Homme_Raw_Hostname} STREQUAL "jaguar")
  SET(Homme_Hostname "Jaguar")
  SET(Homme_Registered_Host TRUE)
ELSE ()
  SET(Homme_Hostname ${Homme_Raw_Hostname})
  SET(Homme_Registered_Host FALSE)
ENDIF()

IF (Homme_Registered_Host)
  MESSAGE(STATUS "Registered Host, reading paths from systemPaths")
  readRegisteredPaths()

ELSE ()
  MESSAGE(STATUS "Host not registered setting hints ")
  determineHintPaths()
  MESSAGE(STATUS "Homme_Hint_Paths=${Homme_Hint_Paths}")
ENDIF ()

# Registered Host specific information

IF (${Homme_Hostname} STREQUAL Yellowstone) 
  SET(Homme_Submission_Type lsf)
ELSEIF (${Homme_Hostname} STREQUAL "Titan")
  SET(Homme_Submission_Type pbs)
ELSE ()
  IF (NOT DEFINED Homme_Submission_Type) 
    SET(Homme_Submission_Type none)
  ENDIF ()
ENDIF () 

IF (NOT ${Homme_Submission_Type} STREQUAL lsf AND
    NOT ${Homme_Submission_Type} STREQUAL pbs AND
    NOT ${Homme_Submission_Type} STREQUAL none)
  #MESSAGE(FATAL_ERROR "Homme_Submission_Type must be one of lsf,pbs,none") 
  MESSAGE(STATUS "Homme_Submission_Type must be one of lsf,pbs,none") 
ENDIF ()  

MESSAGE(STATUS "Homme submission type = ${Homme_Submission_Type}")
