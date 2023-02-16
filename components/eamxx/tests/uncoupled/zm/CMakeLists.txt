include (ScreamUtils)

# Create the test
set (NEED_LIBS zm scream_control scream_share diagnostics)
CreateUnitTest(zm_standalone "zm_standalone.cpp" "${NEED_LIBS}" LABELS "zm;physics;driver"
  MINIMUM_TEST_LEVEL ${SCREAM_TEST_LEVEL_EXPERIMENTAL})

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 24 48)
set (ATM_TIME_STEP 1800)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
