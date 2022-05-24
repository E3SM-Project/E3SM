include (ScreamUtils)

set (NEED_LIBS cld_fraction scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(cld_fraction_standalone "cld_fraction_standalone.cpp" "${NEED_LIBS}" LABELS "cld_fraction;physics;driver")

# Set AD configurable options
set (NUM_STEPS 1)
set (ATM_TIME_STEP 1800)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
