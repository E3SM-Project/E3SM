INCLUDE (ScreamUtils)

SET (zmLib zm)
SET (NEED_LIBS ${zmLib} scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(zm_stand_alone "zm_stand_alone.cpp" "${NEED_LIBS}" LABELS "zm;physics")

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

### Copy zm lookup tables to local data directory
#FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
#CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4
#               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
#CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4
#               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
#CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_universal_constants.inp
#               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
