###########################
# Internal variables      #
###########################

set(OMEGA_PROJECT_NAME            "OmegaOceanModel")
set(OMEGA_EXE_NAME                "omega.exe")
set(OMEGA_LIB_NAME                "OmegaLib")

set(OMEGA_BUILD_MODES             "E3SM" "STANDALONE" "NOT_DEFINED")
set(OMEGA_BUILD_MODE              NOT_DEFINED CACHE STRING "Omega build mode")
set_property(CACHE OMEGA_BUILD_MODE PROPERTY STRINGS ${OMEGA_BUILD_MODES})

set(OMEGA_SOURCE_DIR              ${CMAKE_CURRENT_LIST_DIR})
set(OMEGA_DEFAULT_BUILD_TYPE      Release) # Debug or Release


###########################
# Public variables        #
###########################
macro(setup_common_variables)

  option(OMEGA_DEBUG "Turn on error message throwing (default OFF)." OFF)

  if(NOT DEFINED OMEGA_CXX_FLAGS )
    # initialize cxx flags as an empty list
    set(OMEGA_CXX_FLAGS "")
  endif()

endmacro()

###########################
# Preset Standalone build #
###########################
macro(preset)

  # set CMAKE_CXX_COMPILER from OMEGA_CXX_COMPILER
  if(OMEGA_CXX_COMPILER)
    execute_process(COMMAND which ${OMEGA_CXX_COMPILER}
                    OUTPUT_VARIABLE _OMEGA_CXX_COMPILER
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(CMAKE_CXX_COMPILER ${_OMEGA_CXX_COMPILER})
  endif()

endmacro()

macro(setup_standalone_build)

  setup_common_variables()

  if(NOT DEFINED OMEGA_BUILD_TYPE)
    set(OMEGA_BUILD_TYPE ${OMEGA_DEFAULT_BUILD_TYPE})
  endif()

  if( EXISTS ${OMEGA_SOURCE_DIR}/../../components AND
      EXISTS ${OMEGA_SOURCE_DIR}/../../cime AND
      EXISTS ${OMEGA_SOURCE_DIR}/../../cime_config AND
      EXISTS ${OMEGA_SOURCE_DIR}/../../externals)

    set(E3SM_SOURCE_DIR ${OMEGA_SOURCE_DIR}/../../components)
    set(E3SM_CIME_ROOT ${OMEGA_SOURCE_DIR}/../../cime)
    set(E3SM_CIMECONFIG_ROOT ${OMEGA_SOURCE_DIR}/../../cime_config)
    set(E3SM_EXTERNALS_ROOT ${OMEGA_SOURCE_DIR}/../../externals)

  else()
    # so far, we assume that Omega exists inside of E3SM.
    # However, we leave this else part for later usage.

  endif()

  set(OMEGA_BUILD_MODE "STANDALONE")
  set(OMEGA_BUILD_EXECUTABLE ON)

endmacro()

macro(setup_e3sm_build)

  setup_common_variables()

  set(OMEGA_BUILD_TYPE ${E3SM_DEFAULT_BUILD_TYPE})
  set(E3SM_CIME_ROOT ${CIMEROOT})
  set(E3SM_CIMECONFIG_ROOT ${E3SM_SOURCE_DIR}/../cime_config)
  set(E3SM_EXTERNALS_ROOT ${E3SM_SOURCE_DIR}/../externals)


  #TODO: set OMEGA_ARCH according to E3SM variables
  set(OMEGA_ARCH "NOT_DEFINED")
  set(OMEGA_BUILD_MODE "E3SM")

endmacro()


################################
# Set cmake and YAKL variables #
################################
macro(update_variables)


  if(OMEGA_DEBUG)
    add_compile_definitions(OMEGA_DEBUG)
    add_compile_definitions(LOG_UNBUFFERED_LOGGING="1")
  endif()

  # Set the build type
  set(CMAKE_BUILD_TYPE ${OMEGA_BUILD_TYPE})

  if(OMEGA_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${OMEGA_INSTALL_PREFIX})
  endif()

  if(DEFINED OMEGA_ARCH)

    if(OMEGA_ARCH STREQUAL "NOT_DEFINED")
      set(YAKL_ARCH "")

    else()
      set(YAKL_ARCH ${OMEGA_ARCH})

      if(OMEGA_${YAKL_ARCH}_FLAGS)
        set(YAKL_${YAKL_ARCH}_FLAGS ${OMEGA_${YAKL_ARCH}_FLAGS})
      endif()

    endif()

  else()
    set(YAKL_ARCH "")

  endif()

#  # prints generates all cmake variables
#  get_cmake_property(_variableNames VARIABLES)
#  list (SORT _variableNames)
#  foreach (_variableName ${_variableNames})
#      message(STATUS "${_variableName}=${${_variableName}}")
#  endforeach()

endmacro()



################################
# Verify variable integrity    #
################################
macro(check_setup)

  #message("OMEGA_BUILD_MODE = ${OMEGA_BUILD_MODE}")

  if(OMEGA_BUILD_MODE STREQUAL "E3SM")
    message("*** Omega E3SM-component Build ***")

  elseif(${OMEGA_BUILD_MODE} STREQUAL "STANDALONE")
    message("*** Omega Standalone Build ***")

  else()

    message(FATAL_ERROR "OMEGA_BUILD_MODE is neither E3SM nor STANDALONE.")

  endif()

  if (NOT DEFINED YAKL_ARCH)
    message(FATAL_ERROR "YAKL_ARCH is not defined.")
  endif()

endmacro()


################################
# Prepare output               #
################################
macro(wrap_outputs)

  if(OMEGA_INSTALL_PREFIX)

    install(TARGETS ${OMEGA_LIB_NAME} LIBRARY DESTINATION "${OMEGA_INSTALL_PREFIX}/lib")

    if(OMEGA_BUILD_EXECUTABLE)
      install(TARGETS ${OMEGA_EXE_NAME} RUNTIME DESTINATION "${OMEGA_INSTALL_PREFIX}/bin")
    endif()

  endif()

endmacro()
