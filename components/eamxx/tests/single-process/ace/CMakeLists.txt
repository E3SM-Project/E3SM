include (ScreamUtils)

# Make ACE an optional component that can be disabled
option(SCREAM_ENABLE_ACE "Enable the ACE component" ON)

# Initialize PyTorch found status to FALSE
set(SCREAM_FOUND_PYTORCH FALSE CACHE BOOL "Whether PyTorch was found" FORCE)

# Only proceed if ACE is enabled
if(SCREAM_ENABLE_ACE)
  # Use a configurable path for PyTorch with a default value
  if(NOT DEFINED TORCH_CMAKE_PATH)
    set(TORCH_CMAKE_PATH "/global/cfs/cdirs/e3sm/software/libtorch" 
        CACHE PATH "Path to PyTorch CMake files")
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
else()
  message(STATUS "ACE component is disabled")
endif()
