# Function for handling nl and st gen
function(handle_st_nl_gen NL_GEN_ARGS ST_GEN_ARGS CORE_INPUT_DIR_ARG CORE_BLDDIR_ARG)
  foreach(NL_GEN_ARG IN LISTS NL_GEN_ARGS)
    separate_arguments(ITEMS UNIX_COMMAND "${NL_GEN_ARG}")
    list(GET ITEMS 0 ITEM)
    list(APPEND INPUTS ${ITEM})
    add_custom_command(
      OUTPUT ${CORE_INPUT_DIR_ARG}/${ITEM}
      COMMAND ${CMAKE_BINARY_DIR}/tools/namelist_gen ${CORE_BLDDIR_ARG}/Registry_processed.xml ${NL_GEN_ARG}
      DEPENDS namelist_gen ${CORE_BLDDIR_ARG}/Registry_processed.xml
      WORKING_DIRECTORY ${CORE_INPUT_DIR_ARG}
    )
  endforeach()

  foreach(ST_GEN_ARG IN LISTS ST_GEN_ARGS)
    separate_arguments(ITEMS UNIX_COMMAND "${ST_GEN_ARG}")
    list(GET ITEMS 0 ITEM)
    list(APPEND INPUTS ${ITEM})
    add_custom_command(
      OUTPUT ${CORE_INPUT_DIR_ARG}/${ITEM}
      COMMAND ${CMAKE_BINARY_DIR}/tools/streams_gen ${CORE_BLDDIR_ARG}/Registry_processed.xml ${ST_GEN_ARG}
      DEPENDS streams_gen ${CORE_BLDDIR_ARG}/Registry_processed.xml
      WORKING_DIRECTORY ${CORE_INPUT_DIR_ARG}
    )
  endforeach()

  foreach(INPUT IN LISTS INPUTS)
    add_custom_command(
      OUTPUT ${CORE_BLDDIR_ARG}/${INPUT}
      COMMAND ${CMAKE_COMMAND} -E copy ${CORE_INPUT_DIR_ARG}/${INPUT} ${CORE_BLDDIR_ARG}/${INPUT}
      DEPENDS ${CORE_INPUT_DIR_ARG}/${INPUT}
      WORKING_DIRECTORY ${CORE_BLDDIR_ARG}
    )
  endforeach()
endfunction()

# Function for generating f90 file targets, will add to parent's SOURCES var
function(genf90_targets RAW_SOURCES_ARG INCLUDES_ARG CPPDEFS_ARG NO_PREPROCESS_ARG CORE_INC_DIR_ARG)
  # Add -I to includes so that they can used for cpp command
  foreach(ITEM IN LISTS INCLUDES_ARG)
    list(APPEND INCLUDES_I "-I${ITEM}")
  endforeach()

  # Run all .F files through cpp to generate the f90 file
  foreach(RAW_SOURCE_FILE IN LISTS RAW_SOURCES_ARG)
    get_filename_component(SOURCE_EXT ${RAW_SOURCE_FILE} EXT)
    if ( (SOURCE_EXT STREQUAL ".F" OR SOURCE_EXT STREQUAL ".F90") AND NOT RAW_SOURCE_FILE IN_LIST NO_PREPROCESS_ARG)
      string(REPLACE "${SOURCE_EXT}" ".f90" SOURCE_F90 ${RAW_SOURCE_FILE})
      get_filename_component(DIR_RELATIVE ${SOURCE_F90} DIRECTORY)
      set(DIR_ABSOLUTE ${CMAKE_BINARY_DIR}/${DIR_RELATIVE})
      if (NOT EXISTS ${DIR_ABSOLUTE})
        file(MAKE_DIRECTORY ${DIR_ABSOLUTE})
      endif()
      if (CORE_INC_DIR_ARG)
        add_custom_command (
          OUTPUT ${CMAKE_BINARY_DIR}/${SOURCE_F90}
          COMMAND cpp -P -traditional ${CPPDEFS_ARG} ${INCLUDES_I} -Uvector
          ${CMAKE_CURRENT_SOURCE_DIR}/${RAW_SOURCE_FILE} > ${CMAKE_BINARY_DIR}/${SOURCE_F90}
          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${RAW_SOURCE_FILE} ${CORE_INC_DIR_ARG}/core_variables.inc)
      else()
        add_custom_command (
          OUTPUT ${CMAKE_BINARY_DIR}/${SOURCE_F90}
          COMMAND cpp -P -traditional ${CPPDEFS_ARG} ${INCLUDES_I} -Uvector
          ${CMAKE_CURRENT_SOURCE_DIR}/${RAW_SOURCE_FILE} > ${CMAKE_BINARY_DIR}/${SOURCE_F90}
          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${RAW_SOURCE_FILE})
      endif()
      list(APPEND LOCAL_SOURCES ${CMAKE_BINARY_DIR}/${SOURCE_F90})
    else()
      list(APPEND LOCAL_SOURCES ${RAW_SOURCE_FILE})
    endif()
  endforeach()

  set(SOURCES ${LOCAL_SOURCES} PARENT_SCOPE)

endfunction(genf90_targets)
