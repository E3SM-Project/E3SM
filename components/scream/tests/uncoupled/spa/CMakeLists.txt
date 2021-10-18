include (ScreamUtils)

# Create the test
set (NEED_LIBS spa scream_control scream_share)
CreateUnitTest(spa_stand_alone "spa_stand_alone.cpp" "${NEED_LIBS}" LABELS "spa;physics")

# Set AD configurable options
SetVarDependingOnTestProfile(NUM_STEPS 2 24 48)
set (ATM_TIME_STEP 864000)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)

# Copy input files to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/spa_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/spa_init_ne2np4.nc COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/spa_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/spa_output.yaml COPYONLY)
