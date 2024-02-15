# Gather all sources for a component
#
# Cmake does not have a VPATH concept, so we need relative/absolute paths to source files
# Note: Using absolute paths seems to wreck CMake's ability to do a dep analysis on fortran
# sources and compile them in the right order.
#
# One additional subtley is that, when mkSrcfiles found multiple files with the same basename,
# only the one found first gets compiled. We need to duplicate that behavior here.
function(gather_sources FILEPATH_DIRS_ARG CIMEROOT_ARG)
  set(BASENAME_SET)
  set(SOURCES_RESULT)
  set(GEN_F90_SOURCES_RESULT)

  foreach(DIRSEARCH ${FILEPATH_DIRS_ARG})
    file(GLOB MATCHES RELATIVE "${PROJECT_SOURCE_DIR}" "${DIRSEARCH}/*.[Ffc]" "${DIRSEARCH}/*.[Ff]90" "${DIRSEARCH}/*.cpp" "${DIRSEARCH}/*.F90.in")
    if (MATCHES)
      foreach (MATCH IN LISTS MATCHES)
        get_filename_component(BASENAME ${MATCH} NAME)
        list(FIND BASENAME_SET ${BASENAME} BASENAME_WAS_FOUND)
        if (BASENAME_WAS_FOUND EQUAL -1)
          list(APPEND SOURCES_RESULT ${MATCH})
          list(APPEND BASENAME_SET ${BASENAME})
        else()
          message("Warning: Skipping repeated base filename ${BASENAME} for ${MATCH}")
        endif()
      endforeach()
    endif()
  endforeach()

  foreach(SOURCE_FILE IN LISTS SOURCES_RESULT)
    get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)
    if (SOURCE_EXT STREQUAL ".F90.in")
      string(REPLACE ".in" "" SOURCE_NO_IN ${SOURCE_FILE})
      list(APPEND GEN_F90_SOURCES_RESULT ${SOURCE_NO_IN})
      list(APPEND SOURCES_RESULT ${SOURCE_NO_IN})
      list(REMOVE_ITEM SOURCES_RESULT ${SOURCE_FILE})
    endif()
  endforeach()

  # Return data to parent
  set(SOURCES_RESULT ${SOURCES_RESULT} PARENT_SCOPE)
  set(GEN_F90_SOURCES_RESULT ${GEN_F90_SOURCES_RESULT} PARENT_SCOPE)

endfunction()

function(e3sm_get_source_property FILE_ARG PROPERTY_NAME)
  if (FILE_ARG MATCHES "${CMAKE_BINARY_DIR}/.*") # is generated
    set(REAL_FILE ${FILE_ARG})
  else()
    set(REAL_FILE "${PROJECT_SOURCE_DIR}/${FILE_ARG}")
  endif()
  get_property(RESULT SOURCE ${REAL_FILE} PROPERTY ${PROPERTY_NAME})

  # Return data to parent
  set(PROPERTY_RESULT ${RESULT} PARENT_SCOPE)
endfunction()

function(e3sm_set_source_property FILE_ARG PROPERTY_NAME PROPERTY_VALUE IS_APPEND)
  if (FILE_ARG MATCHES "${CMAKE_BINARY_DIR}/.*") # is generated
    set(REAL_FILE ${FILE_ARG})
  else()
    if (NOT EXISTS ${PROJECT_SOURCE_DIR}/${FILE_ARG})
      message(FATAL_ERROR "Trying to set property on non-existent source: ${FILE_ARG}, looked for ${PROJECT_SOURCE_DIR}/${FILE_ARG}")
    endif()
    set(REAL_FILE "${PROJECT_SOURCE_DIR}/${FILE_ARG}")
  endif()
  if (IS_APPEND)
    set_property(SOURCE ${REAL_FILE} APPEND_STRING PROPERTY ${PROPERTY_NAME} "${PROPERTY_VALUE}")
  else()
    set_property(SOURCE ${REAL_FILE}               PROPERTY ${PROPERTY_NAME} "${PROPERTY_VALUE}")
  endif()
endfunction()

# Add compile flags for a file. Expects a filepath relative to E3SM/components if it
# is not a generated file.
function(e3sm_add_flags FILE_ARG FLAGS_ARG)
  e3sm_set_source_property(${FILE_ARG} COMPILE_FLAGS " ${FLAGS_ARG}" TRUE)
endfunction()

function(e3sm_deoptimize_file FILE_ARG)
  e3sm_get_source_property(${FILE_ARG} LANGUAGE)
  # Adding all DEBUG flags may be overkill
  #e3sm_add_flags(${FILE_ARG} "${CMAKE_${PROPERTY_RESULT}_FLAGS_DEBUG}")
  e3sm_add_flags(${FILE_ARG} "-O0")
  if (COMPILER STREQUAL "nvidia")
    e3sm_add_flags(${FILE_ARG} "-Mnofma")
  endif()
endfunction()
