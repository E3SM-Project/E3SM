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
    file(GLOB MATCHES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/../.." "${DIRSEARCH}/*.[Ffc]" "${DIRSEARCH}/*.[Ff]90" "${DIRSEARCH}/*.cpp" "${DIRSEARCH}/*.F90.in")
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

# Add compile flags for a file. Expects a filepath relative to E3SM/components if it
# is not a generated file.
function(e3sm_add_flags FILE_ARG FLAGS_ARG)
  if (FILE_ARG MATCHES "${CMAKE_BINARY_DIR}/.*") # is generated
    set(REAL_FILE ${FILE_ARG})
  else()
    if (NOT EXISTS ${PROJECT_SOURCE_DIR}/${FILE_ARG})
      message(FATAL_ERROR "Trying to set flags on non-existent source: ${FILE_ARG}, looked for ${PROJECT_SOURCE_DIR}/${FILE_ARG}")
    endif()
    set(REAL_FILE "${SOURCE_PATH}/${FILE_ARG}")
  endif()
  set_property(SOURCE ${REAL_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FLAGS_ARG} ")
endfunction()

# Remove compile flags for a file. Expects a filepath relative to E3SM/components if it
# is not a generated file.
function(e3sm_remove_flags FILE_ARG FLAGS_ARG)
  if (FILE_ARG MATCHES "${CMAKE_BINARY_DIR}/.*") # is generated
    set(REAL_FILE ${FILE_ARG})
  else()
    if (NOT EXISTS ${PROJECT_SOURCE_DIR}/${FILE_ARG})
      message(FATAL_ERROR "Trying to set flags on non-existent source: ${FILE_ARG}")
    endif()
    set(REAL_FILE "${SOURCE_PATH}/${FILE_ARG}")
  endif()

  get_property(ITEM_FLAGS SOURCE ${REAL_FILE} PROPERTY COMPILE_FLAGS)
  string(REPLACE "${FLAGS_ARG}" "" ITEM_FLAGS "${ITEM_FLAGS}")
  set_property(SOURCE ${REAL_FILE} PROPERTY COMPILE_FLAGS "${ITEM_FLAGS}")
endfunction()

function(e3sm_deoptimize_file FILE_ARG FFLAGS_NOOPT)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    # PGI does not support bulk-disabling of optimization by appending -O0,
    # we have to remove the optimization flags first

    # Until we know which particular flags are related to optimization, we have to guess.
    e3sm_remove_flags(${FILE_ARG} "-O1")
    e3sm_remove_flags(${FILE_ARG} "-O2")
    e3sm_remove_flags(${FILE_ARG} "-O3")
  endif()
  e3sm_add_flags(${FILE_ARG} "${FFLAGS_NOOPT}")
endfunction()
