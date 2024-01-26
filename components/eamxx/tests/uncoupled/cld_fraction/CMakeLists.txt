include (ScreamUtils)

# Test atmosphere processes
CreateADUnitTest(cld_fraction_standalone
  LIBS cld_fraction
  LABELS cld_fraction physics)

# Set AD configurable options
set (NUM_STEPS 1)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
