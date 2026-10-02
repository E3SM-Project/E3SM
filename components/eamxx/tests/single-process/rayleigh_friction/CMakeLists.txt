include (ScreamUtils)

# Set common AD configurable options
set (NUM_STEPS 2)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Create executable for RayleighFriction standalone
CreateADUnitTestExec(rayleigh_friction_standalone
  LIBS rayleigh_friction)

# Configure yaml files to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

CreateUnitTestFromExec(rayleigh_friction_standalone rayleigh_friction_standalone
  EXE_ARGS "--args -ifile=input.yaml"
  LABELS rayleigh_friction physics)
