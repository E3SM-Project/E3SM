# For now, scream will rely on it's own build of kokkos rather than the
# one in sharedlib.

function(build_scream)

 string(FIND "${CAM_CONFIG_OPTS}" "-crm samxx" HAS_SAMXX)
 if (NOT HAS_SAMXX EQUAL -1)
   # The following is for the SAMXX code:
   set(USE_SAMXX TRUE CACHE INTERNAL "use samxx")
 endif()

 if (USE_SAMXX  OR COMP_NAMES MATCHES ".*scream.*") 
    # Include machine file here
    message("Found scream component")
    if (USE_CUDA)
       include(${CMAKE_SOURCE_DIR}/eam/src/physics/crm/scream/cmake/machine-files/${MACH}.cmake)
    else()
       include(${CMAKE_SOURCE_DIR}/eam/src/physics/crm/scream/cmake/machine-files/${MACH}_cpu.cmake)
    endif()

    set(SCREAM_SOURCE_DIR ${CMAKE_SOURCE_DIR}/eam/src/physics/crm/scream)
    add_subdirectory(${SCREAM_SOURCE_DIR} ./scream)
    set(SCREAM_SOURCE_DIR ${CMAKE_SOURCE_DIR}/eam/src/physics/crm/scream/src CACHE STRING "Scream Source Dir")
    set(SCREAM_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/eam/src/physics/crm/scream CACHE STRING "Scream Binary Dir")
 endif()

endfunction(build_scream)
