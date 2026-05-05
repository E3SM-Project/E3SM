# Emulator Components build configuration
# Handles all emulator components: EMULATORATM, EMULATOROCN
# Similar pattern to build_eamxx.cmake

# Capture this file's directory at include-time.
# CMAKE_CURRENT_LIST_DIR inside a function body is evaluated at *call* time
# (i.e. it becomes the caller's directory), so we must snapshot it here.
set(_EMULATOR_BUILD_COMPS_DIR ${CMAKE_CURRENT_LIST_DIR})

function(build_emulator_comps)

  # Check if any emulator component is being built
  set(EMULATORATM_FOUND FALSE)
  set(EMULATOROCN_FOUND FALSE)

  if (COMP_NAMES MATCHES ".*emulatoratm.*")
    set(EMULATORATM_FOUND TRUE)
  endif()
  if (COMP_NAMES MATCHES ".*emulatorocn.*")
    set(EMULATOROCN_FOUND TRUE)
  endif()


  # Only proceed if at least one emulator component is used
  if (EMULATORATM_FOUND OR EMULATOROCN_FOUND)
    
    message(STATUS "")
    message(STATUS "=================================================================")
    message(STATUS "  Building Emulator Components Framework")
    message(STATUS "=================================================================")
    message(STATUS "  Components: EMULATORATM=${EMULATORATM_FOUND} EMULATOROCN=${EMULATOROCN_FOUND}")

    include(${_EMULATOR_BUILD_COMPS_DIR}/common_setup.cmake)

    #---------------------------------------------------------------------------
    # Build emulator_comps using add_subdirectory
    #---------------------------------------------------------------------------
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

    #---------------------------------------------------------------------------
    # Export variables for build_model.cmake to link against emulator libs
    #---------------------------------------------------------------------------
    # These variables will be used in build_model.cmake when building the 
    # emulator components to link against the emulator libraries
    set(EMULATOR_COMPS_BUILT TRUE CACHE INTERNAL "Emulator components were built")
    set(EMULATOR_COMMON_LIB emulator_common CACHE INTERNAL "Common emulator library")
    set(EMULATOR_DRIVER_LIB emulator_driver CACHE INTERNAL "Emulator driver library")
    set(EMULATORATM_LIB emulatoratm CACHE INTERNAL "EMULATORATM library")
    set(EMULATOROCN_LIB emulatorocn CACHE INTERNAL "EMULATOROCN library")

    # Alias emulatoratm as the 'atm' component target so that build_model.cmake
    # links against libemulatoratm.a, which in the full E3SM build includes
    # emulator_create (atm_factory.cpp) directly — no circular dep on
    # emulator_driver.
    if(EMULATORATM_FOUND AND TARGET emulatoratm)
      add_library(atm ALIAS emulatoratm)
    endif()

    #---------------------------------------------------------------------------
    # Export list of emulator components (for SKIP_COMPS in CMakeLists.txt)
    #---------------------------------------------------------------------------
    set(_EMULATOR_COMP_NAMES "")
    if(EMULATORATM_FOUND)
      list(APPEND _EMULATOR_COMP_NAMES "emulatoratm")
    endif()
    if(EMULATOROCN_FOUND)
      list(APPEND _EMULATOR_COMP_NAMES "emulatorocn")
    endif()

    set(EMULATOR_COMP_NAMES ${_EMULATOR_COMP_NAMES} PARENT_SCOPE)
    set(EMULATOR_COMP_NAMES ${_EMULATOR_COMP_NAMES} CACHE INTERNAL "List of emulator component names")

    message(STATUS "  Emulator components to skip in build_model: ${_EMULATOR_COMP_NAMES}")

    message(STATUS "=================================================================")
    message(STATUS "  Emulator Components configuration complete")
    message(STATUS "=================================================================")
    message(STATUS "")
  endif()

endfunction(build_emulator_comps)
