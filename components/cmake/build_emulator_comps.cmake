# Emulator Components build configuration
# Handles the emulator framework and component shims.
#
# Important design point:
#   CMake should only decide whether the emulator framework is built.
#   It should not decide which emulator implementation is active.
#   Runtime configuration and the emulator factory decide that.

# Capture this file's directory at include-time.
# CMAKE_CURRENT_LIST_DIR inside a function body is evaluated at call time.
set(_EMULATOR_BUILD_COMPS_DIR ${CMAKE_CURRENT_LIST_DIR})

function(build_emulator_comps)

  # Detect whether any emulator component was requested.
  set(EMULATOR_ENABLED FALSE)
  set(_EMULATOR_COMP_NAMES "")

  foreach(_comp IN LISTS COMP_NAMES)
    if(_comp MATCHES "^emulator")
      set(EMULATOR_ENABLED TRUE)
      list(APPEND _EMULATOR_COMP_NAMES "${_comp}")
    endif()
  endforeach()

  if(NOT EMULATOR_ENABLED)
    return()
  endif()

  message(STATUS "")
  message(STATUS "=================================================================")
  message(STATUS "  Building Emulator Components Framework")
  message(STATUS "=================================================================")
  message(STATUS "  Emulator component requests: ${_EMULATOR_COMP_NAMES}")

  include(${_EMULATOR_BUILD_COMPS_DIR}/common_setup.cmake)

  # ---------------------------------------------------------------------------
  # Build the complete emulator framework.
  # ---------------------------------------------------------------------------
  set(EMULATOR_COMPS_DIR ${_EMULATOR_BUILD_COMPS_DIR}/../emulators)

  if(EXISTS ${EMULATOR_COMPS_DIR}/CMakeLists.txt)
    message(STATUS "  Including emulator components from:")
    message(STATUS "    ${EMULATOR_COMPS_DIR}")
    add_subdirectory(${EMULATOR_COMPS_DIR} ${CMAKE_BINARY_DIR}/emulators)
  else()
    message(FATAL_ERROR
      "Emulator components directory not found!\n"
      "  Expected: ${EMULATOR_COMPS_DIR}/CMakeLists.txt\n"
      "  Current source dir: ${CMAKE_SOURCE_DIR}\n"
      "Please ensure the emulator components are present in the source tree.")
  endif()

  # ---------------------------------------------------------------------------
  # Export variables for build_model.cmake.
  #
  # These are framework/build integration variables only. Runtime configuration
  # chooses the actual emulator behavior through the emulator factory.
  # ---------------------------------------------------------------------------
  set(EMULATOR_COMPS_BUILT TRUE CACHE INTERNAL "Emulator components were built")
  set(EMULATOR_COMMON_LIB emulator_common CACHE INTERNAL "Common emulator library")
  set(EMULATOR_DRIVER_LIB emulator_driver CACHE INTERNAL "Emulator driver library")

  if(TARGET emulatoratm)
    set(EMULATORATM_LIB emulatoratm CACHE INTERNAL "EMULATORATM component shim library")
  endif()

  if(TARGET emulatorocn)
    set(EMULATOROCN_LIB emulatorocn CACHE INTERNAL "EMULATOROCN component shim library")
  endif()

  # ---------------------------------------------------------------------------
  # Component aliases.
  #
  # These aliases satisfy E3SM's expected component target names. They should be
  # thin integration shims; emulator behavior itself should come from runtime
  # configuration and the emulator factory.
  # ---------------------------------------------------------------------------
  if("emulatoratm" IN_LIST _EMULATOR_COMP_NAMES AND TARGET emulatoratm)
    add_library(atm ALIAS emulatoratm)
  endif()

  if("emulatorocn" IN_LIST _EMULATOR_COMP_NAMES AND TARGET emulatorocn)
    add_library(ocn ALIAS emulatorocn)
  endif()

  # ---------------------------------------------------------------------------
  # Export emulator component names so the normal component build can skip them.
  # ---------------------------------------------------------------------------
  set(EMULATOR_COMP_NAMES ${_EMULATOR_COMP_NAMES} PARENT_SCOPE)
  set(EMULATOR_COMP_NAMES ${_EMULATOR_COMP_NAMES} CACHE INTERNAL "List of emulator component names")

  message(STATUS "  Emulator components to skip in build_model: ${_EMULATOR_COMP_NAMES}")
  message(STATUS "=================================================================")
  message(STATUS "  Emulator Components configuration complete")
  message(STATUS "=================================================================")
  message(STATUS "")

endfunction()
