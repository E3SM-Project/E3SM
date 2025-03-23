include (ScreamUtils)

# Find PyTorch and set up dependencies
set(CMAKE_PREFIX_PATH "/global/u1/m/mahf708/.micromamba/envs/fme/lib/python3.10/site-packages/torch/share/cmake/Torch")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
find_package(Torch REQUIRED)

# Test atmosphere processes: added physics_share, scream_share, and ekat libraries
CreateADUnitTest(ace_standalone
  LIBS ace physics_share scream_share scream ${TORCH_LIBRARIES} ekat
  COMPILER_DEFS EAMXX_HAS_ACE
  INCLUDE_DIRS ${TORCH_INCLUDE_DIRS}
  LINKER_FLAGS "-Wl,--copy-dt-needed-entries"
  LABELS ace physics)

# Set AD configurable options
set (NUM_STEPS 1)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
