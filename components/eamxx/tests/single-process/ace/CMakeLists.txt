include (ScreamUtils)

# Only build the ACE test if PyTorch was found and ACE is enabled
if(SCREAM_ENABLE_ACE AND SCREAM_FOUND_PYTORCH)
  CreateADUnitTest(ace_standalone
    LABELS ace physics
    LIBS ace physics_share scream_share ekat
  )
  
  # Set AD configurable options
  set (NUM_STEPS 1)
  set (ATM_TIME_STEP 1800)
  set (RUN_T0 2021-10-12-45000)
  
  # Configure yaml input file to run directory
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
else()
  message(STATUS "Skipping ACE standalone test since ACE is not enabled or PyTorch was not found")
endif()
