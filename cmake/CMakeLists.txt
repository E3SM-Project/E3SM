# This file is used primarily to provide macro definitions and set the default locations for dependencies on supported systems


# Determine the host name
EXECUTE_PROCESS(COMMAND uname -n
  RESULT_VARIABLE Homme_result
  OUTPUT_VARIABLE Homme_output
  ERROR_VARIABLE Homme_error
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

IF (Homme_result EQUAL 0 AND Homme_error STREQUAL "")
  SET(Homme_Hostname ${Homme_output})
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

MESSAGE(STATUS "Homme hostname = ${Homme_Hostname}")
MESSAGE(STATUS "Homme OS = ${Homme_OS}")

# Now parse the host name to see if it matches yslogin*
IF (${Homme_Hostname} STREQUAL "yslogin3")
  MESSAGE(STATUS "Hostname is yslogin3")
ENDIF()

IF (${Homme_Hostname} MATCHES "yslogin.")
  MESSAGE(STATUS "Hostname is yslogin")
ENDIF()
