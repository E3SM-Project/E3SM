INCLUDE (ScreamUtils)

SET (shocLib shoc)
SET (NEED_LIBS ${shocLib} scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(shoc_stand_alone "shoc_stand_alone.cpp" "${NEED_LIBS}" LABELS "shoc;physics")

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/shoc_standalone_ne4_ne4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_standalone_ne4_ne4.nc COPYONLY)
