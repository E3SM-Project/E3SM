include(ScreamUtils)

# Set GPU_OR_CPU variable (to use in input.yaml for checkpoint)
if (Kokkos_ENABLE_CUDA)
  set(GPU_OR_CPU "" CACHE STRING "GPU or CPU version of PyTorch.")
else()
  set(GPU_OR_CPU "_cpu" CACHE STRING "GPU or CPU version of PyTorch.")
endif()

# Create the test with the determined torch target
CreateADUnitTest(ace_standalone
  LIBS ace
  LABELS ace physics)

# Set AD configurable options
set (NUM_STEPS 2)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
              ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
# Configure yaml output file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
              ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)
