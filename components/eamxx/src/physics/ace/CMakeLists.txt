set(ACE_SRCS
  eamxx_ace_process_interface.cpp
)

# Make ACE an optional component that can be disabled
option(SCREAM_ENABLE_ACE "Enable the ACE component" ON)

# Initialize PyTorch found status to FALSE
set(SCREAM_FOUND_PYTORCH FALSE CACHE BOOL "Whether PyTorch was found" FORCE)

# Only proceed if ACE is enabled
if(SCREAM_ENABLE_ACE)
  # Use a configurable path for PyTorch with a default value
  if(NOT DEFINED TORCH_CMAKE_PATH)
    set(TORCH_CMAKE_PATH "/global/u1/m/mahf708/.micromamba/envs/fme/lib/python3.10/site-packages/torch/share/cmake/Torch" 
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
    add_library(ace ${ACE_SRCS})
    target_compile_definitions(ace PUBLIC EAMXX_HAS_ACE)
    target_include_directories(ace SYSTEM PUBLIC ${TORCH_INCLUDE_DIRS})
    
    target_link_libraries(ace physics_share scream_share ${TORCH_LIBRARIES})
    
    if (TARGET eamxx_physics)
      target_link_libraries(eamxx_physics INTERFACE ace)
    endif()
  else()
    message(WARNING "PyTorch not found at ${TORCH_CMAKE_PATH}. ACE component will not be built.")
    message(STATUS "To install PyTorch, follow instructions at https://pytorch.org/get-started/locally/")
    message(STATUS "Or specify custom PyTorch location with -DTORCH_CMAKE_PATH=/path/to/torch/cmake")
    
    # Set a variable to indicate that ACE is not available
    set(SCREAM_ENABLE_ACE OFF CACHE BOOL "Enable the ACE component" FORCE)
  endif()
else()
  message(STATUS "ACE component is disabled")
endif()
