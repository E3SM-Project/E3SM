# Emulator Components build configuration
# Handles all emulator components: EATM, EOCN
# Similar pattern to build_eamxx.cmake

function(build_emulator_comps)

  # Check if any emulator component is being built
  set(EATM_FOUND FALSE)
  set(EOCN_FOUND FALSE)
  set(EICE_FOUND FALSE)

  if (COMP_NAMES MATCHES ".*eatm.*")
    set(EATM_FOUND TRUE)
  endif()
  if (COMP_NAMES MATCHES ".*eocn.*")
    set(EOCN_FOUND TRUE)
  endif()


  # Only proceed if at least one emulator component is used
  if (EATM_FOUND OR EOCN_FOUND)
    
    message(STATUS "")
    message(STATUS "=================================================================")
    message(STATUS "  Building Emulator Components Framework")
    message(STATUS "=================================================================")
    message(STATUS "  Components: EATM=${EATM_FOUND} EOCN=${EOCN_FOUND}")

    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

    #---------------------------------------------------------------------------
    # Build emulator_comps using add_subdirectory
    #---------------------------------------------------------------------------
    set(EMULATOR_COMPS_DIR ${CMAKE_SOURCE_DIR}/emulators)
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
    set(EMULATOR_EATM_LIB eatm CACHE INTERNAL "EATM library")
    set(EMULATOR_EOCN_LIB eocn CACHE INTERNAL "EOCN library")

    # Create aliases for common component class names so build_model.cmake works
    if(EATM_FOUND AND TARGET eatm)
      add_library(atm ALIAS eatm)
    endif()
    if(EOCN_FOUND AND TARGET eocn)
      add_library(ocn ALIAS eocn)
    endif()

    #---------------------------------------------------------------------------
    # Export list of emulator components (for SKIP_COMPS in CMakeLists.txt)
    #---------------------------------------------------------------------------
    set(_EMULATOR_COMP_NAMES "")
    if(EATM_FOUND)
      list(APPEND _EMULATOR_COMP_NAMES "eatm")
    endif()
    if(EOCN_FOUND)
      list(APPEND _EMULATOR_COMP_NAMES "eocn")
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
