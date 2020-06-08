macro (EkatConfigFile CONFIG_FILE_IN CONFIG_FILE_C CONFIG_FILE_F90)
  set(options OPTIONAL AT_ONLY)
  set(oneValueArgs)
  set(multiValueArgs)

  cmake_parse_arguments(EKAT_CONFIGURE_FILE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # Generate temporary config file
  if (EKAT_CONFIGURE_FILE_AT_ONLY)
    configure_file (${CONFIG_FILE_IN} ${CONFIG_FILE_C}.tmp @ONLY)
  else()
    configure_file (${CONFIG_FILE_IN} ${CONFIG_FILE_C}.tmp)
  endif()

  # Assume by default that config file is out of date
  set (OUT_OF_DATE TRUE)

  # If config file in binary dir exists, we check whether the new one would be different
  if (EXISTS ${CONFIG_FILE_C})

    # We rely on FILE macro rather than running diff, since it is
    # more portable (guaranteed to work regardless of underlying system)
    file (READ ${CONFIG_FILE_C} CONFIG_FILE_C_STR)
    file (READ ${CONFIG_FILE_C}.tmp CONFIG_FILE_C_TMP_STR)

    if (${CONFIG_FILE_C_STR} STREQUAL ${CONFIG_FILE_C_TMP_STR})
      # config file was present and appears unchanged
      set (OUT_OF_DATE FALSE)
    endif()

    FILE (REMOVE ${CONFIG_FILE_C}.tmp)
  endif ()

  # If out of date (either missing or different), adjust
  if (OUT_OF_DATE)

    # Run the configure macro
    configure_file (${CONFIG_FILE_IN} ${CONFIG_FILE_C})

    # run sed to change '/*...*/' comments into '!/*...*/'
    execute_process(COMMAND sed "s;^/;!/;g"
                    WORKING_DIRECTORY ${EKAT_BINARY_DIR}
                    INPUT_FILE ${CONFIG_FILE_C}
                    OUTPUT_FILE ${CONFIG_FILE_F90})

    # do the same for '//...' comments (turn them into '! ...'
    execute_process(COMMAND sed "s;^//;!;g"
                    WORKING_DIRECTORY ${EKAT_BINARY_DIR}
                    INPUT_FILE ${CONFIG_FILE_C}
                    OUTPUT_FILE ${CONFIG_FILE_F90})
  endif()

endmacro (EkatConfigFile)
