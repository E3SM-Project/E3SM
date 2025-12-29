# Emulator Components build configuration
# Handles all emulator components: EATM, EOCN, EICE (future)
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
  if (COMP_NAMES MATCHES ".*eice.*")
    set(EICE_FOUND TRUE)
  endif()

  # Only proceed if at least one emulator component is used
  if (EATM_FOUND OR EOCN_FOUND OR EICE_FOUND)
    
    message(STATUS "")
    message(STATUS "=================================================================")
    message(STATUS "  Building Emulator Components Framework")
    message(STATUS "=================================================================")
    message(STATUS "  Components: EATM=${EATM_FOUND} EOCN=${EOCN_FOUND} EICE=${EICE_FOUND}")
    
    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

    #---------------------------------------------------------------------------
    # yaml-cpp: Required for YAML configuration parsing
    #---------------------------------------------------------------------------
    message(STATUS "  Looking for yaml-cpp...")
    find_package(yaml-cpp QUIET)
    
    if(yaml-cpp_FOUND)
      message(STATUS "    Found system yaml-cpp: ${yaml-cpp_DIR}")
      # Export the target for use by emulator_comps
      set(EMULATOR_YAML_CPP_TARGET yaml-cpp::yaml-cpp CACHE INTERNAL "yaml-cpp target")
      set(EMULATOR_YAML_CPP_LIBS yaml-cpp CACHE INTERNAL "yaml-cpp library for final linking")
    else()
      # Build yaml-cpp from EKAT submodule
      set(YAML_CPP_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../externals/ekat/extern/yaml-cpp)
      if(EXISTS ${YAML_CPP_SOURCE_DIR}/CMakeLists.txt)
        message(STATUS "    Building yaml-cpp from EKAT submodule:")
        message(STATUS "      Source: ${YAML_CPP_SOURCE_DIR}")
        
        # Set options before adding subdirectory
        set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
        set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
        set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
        set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
        
        set(YAML_CPP_BUILD_DIR ${CMAKE_BINARY_DIR}/externals/yaml-cpp)
        add_subdirectory(${YAML_CPP_SOURCE_DIR} ${YAML_CPP_BUILD_DIR})
        
        # Create alias if not exists
        if(NOT TARGET yaml-cpp::yaml-cpp)
          add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp)
        endif()
        
        set(EMULATOR_YAML_CPP_TARGET yaml-cpp CACHE INTERNAL "yaml-cpp target")
        set(EMULATOR_YAML_CPP_LIBS yaml-cpp CACHE INTERNAL "yaml-cpp library for final linking")
        message(STATUS "      Built yaml-cpp target: ${EMULATOR_YAML_CPP_TARGET}")
      else()
        message(FATAL_ERROR "yaml-cpp not found!\n"
                "  Checked: system (find_package)\n"
                "  Checked: EKAT submodule at ${YAML_CPP_SOURCE_DIR}\n"
                "  Please ensure yaml-cpp is available.")
      endif()
    endif()

    #---------------------------------------------------------------------------
    # Build emulator_comps using add_subdirectory
    #---------------------------------------------------------------------------
    set(EMULATOR_COMPS_DIR ${CMAKE_SOURCE_DIR}/emulator_comps)
    if(EXISTS ${EMULATOR_COMPS_DIR}/CMakeLists.txt)
      message(STATUS "  Including emulator_comps from:")
      message(STATUS "    ${EMULATOR_COMPS_DIR}")
      add_subdirectory(${EMULATOR_COMPS_DIR} ${CMAKE_BINARY_DIR}/emulator_comps)
    else()
      message(FATAL_ERROR "emulator_comps directory not found at ${EMULATOR_COMPS_DIR}")
    endif()

    #---------------------------------------------------------------------------
    # Export variables for build_model.cmake to link against emulator libs
    #---------------------------------------------------------------------------
    # These variables will be used in build_model.cmake when building the 
    # eatm/eocn/eice components to link against the emulator libraries
    set(EMULATOR_COMPS_BUILT TRUE CACHE INTERNAL "Emulator components were built")
    set(EMULATOR_COMMON_LIB emulator_common CACHE INTERNAL "Common emulator library")
    set(EMULATOR_EATM_LIB eatm CACHE INTERNAL "EATM library")
    
    #---------------------------------------------------------------------------
    # Export list of emulator components (for SKIP_COMPS in CMakeLists.txt)
    #---------------------------------------------------------------------------
    set(EMULATOR_COMP_NAMES "eatm" CACHE INTERNAL "List of emulator component names")
    # Future: append eocn, eice as they're added
    if(EOCN_FOUND)
      list(APPEND EMULATOR_COMP_NAMES "eocn")
    endif()
    if(EICE_FOUND)
      list(APPEND EMULATOR_COMP_NAMES "eice")
    endif()
    
    message(STATUS "  Emulator components to skip in build_model: ${EMULATOR_COMP_NAMES}")
    
    message(STATUS "=================================================================")
    message(STATUS "  Emulator Components configuration complete")
    message(STATUS "=================================================================")
    message(STATUS "")
  endif()

endfunction(build_emulator_comps)
