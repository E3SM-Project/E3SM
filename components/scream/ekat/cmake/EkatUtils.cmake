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

macro (EkatDisableAllWarning targetName)

  if (NOT TARGET ${targetName})
    message (FATAL_ERROR "Error! Cannot disable warnings for target ${targetName}; it is not built by this project.")
  endif ()

  # Add flags to ignore warnings to the target, for all Ekat-supported languages (C, CXX, Fortran)
  # Make the flag compiler-dependent. Notice that only one of the $<$<C_COMPILER_ID:blah>: "blahblah">
  # will expand to anything at all, so this is ok.
  # Note: even if a compiler collection (usually) has the same flag for all languages, we still
  #       add the flag separately for each langauge, since the user MAY be using different compilers
  #       for different langauges (e.g., icpc and gfortran).
  target_compile_options (${targetName} PRIVATE
    $<$<COMPILE_LANGUAGE:C>:$<$<C_COMPILER_ID:GNU>:-w> $<$<C_COMPILER_ID:Intel>: -warn>>)

  target_compile_options (${targetName} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:$<$<Fortran_COMPILER_ID:GNU>:-w> $<$<Fortran_COMPILER_ID:Intel>: -warn>>)
  target_compile_options (${targetName} PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CXX_COMPILER_ID:GNU>:-w> $<$<CXX_COMPILER_ID:Intel>: -warn>>)

endmacro (EkatDisableAllWarning)

