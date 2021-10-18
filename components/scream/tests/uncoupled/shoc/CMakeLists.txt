include (ScreamUtils)

# Create the test
set (NEED_LIBS shoc scream_control scream_share)
CreateUnitTest(shoc_standalone "shoc_standalone.cpp" "${NEED_LIBS}" LABELS "shoc;physics")

# Set AD configurable options
SetVarDependingOnTestProfile(NUM_STEPS 2 24 48)
set (ATM_TIME_STEP 1800)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)

# Copy input files to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shoc_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_init_ne2np4.nc COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shoc_standalone_output.yaml 
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_standalone_output.yaml COPYONLY)
