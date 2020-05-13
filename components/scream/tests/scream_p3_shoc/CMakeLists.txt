INCLUDE (ScreamUtils)

SET (p3Lib p3)
SET (shocLib shoc)
SET (NEED_LIBS ${p3Lib} ${shocLib} scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(p3_shoc_coupled "p3_shoc_coupled.cpp" "${NEED_LIBS}" LABELS "p3;shoc;physics")

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

# Copy p3 lookup tables to local data directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_universal_constants.inp
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
