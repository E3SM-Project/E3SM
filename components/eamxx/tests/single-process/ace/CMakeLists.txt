# Initialize PyTorch found status to FALSE
set(SCREAM_FOUND_PYTORCH FALSE CACHE BOOL "Whether PyTorch was found" FORCE)

# Use a configurable path for PyTorch with a default value
if(NOT DEFINED TORCH_CMAKE_PATH)
  if (Kokkos_ENABLE_CUDA)
    set(GPU_OR_CPU "" CACHE STRING "GPU or CPU version of PyTorch.")
    set(TORCH_CMAKE_PATH "/global/cfs/cdirs/e3sm/software/libtorch/libtorch-cxx11-abi-shared-with-deps-2.6.0+cu124" 
        CACHE PATH "Path to PyTorch CMake files")
    message(STATUS "Will try to use GPU version of PyTorch at ${TORCH_CMAKE_PATH} with ${GPU_OR_CPU}")
  else()
    set(GPU_OR_CPU "_cpu" CACHE STRING "GPU or CPU version of PyTorch.")
    set(TORCH_CMAKE_PATH "/global/cfs/cdirs/e3sm/software/libtorch/libtorch-cxx11-abi-shared-with-deps-2.6.0+cpu" 
        CACHE PATH "Path to PyTorch CMake files")
    message(STATUS "Will try to use CPU version of PyTorch at ${TORCH_CMAKE_PATH} with ${GPU_OR_CPU}")
  endif()
  # TODO add support for ROCm
endif()

# Try to find PyTorch, but don't require it
set(CMAKE_PREFIX_PATH ${TORCH_CMAKE_PATH})
find_package(Torch QUIET)

# Only build ACE if PyTorch is found
if(Torch_FOUND)
  message(STATUS "Found PyTorch: Building ACE component")
  # Set variable to indicate PyTorch was found, make it visible to parent scope
  set(SCREAM_FOUND_PYTORCH TRUE CACHE BOOL "Whether PyTorch was found" FORCE)
  
  include(ScreamUtils)
  # Add the ACE test
  # Only build the ACE test if PyTorch was found and ACE is enabled

  CreateADUnitTest(ace_standalone
    LIBS ${TORCH_LIBRARIES} ace
    LABELS ace physics)
  
  target_link_libraries(ace_standalone PUBLIC
    ace
    ${TORCH_LIBRARIES}
  )

  if (Kokkos_ENABLE_CUDA)
    set(GPU_OR_CPU "" CACHE STRING "GPU or CPU version of PyTorch.")
  else()
    set(GPU_OR_CPU "_cpu" CACHE STRING "GPU or CPU version of PyTorch.")
  endif()

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

else()
  message(STATUS "Skipping ACE standalone test since ACE is not enabled or PyTorch was not found")
endif()
