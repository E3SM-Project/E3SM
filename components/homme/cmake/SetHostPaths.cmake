
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
ELSEIF (${Homme_Raw_Hostname} MATCHES "titan-ext." OR ${Homme_Raw_Hostname} MATCHES "titan-login.")
  SET(Homme_Hostname "Titan")
ELSE ()
  SET(Homme_Hostname ${Homme_Raw_Hostname})
ENDIF()

