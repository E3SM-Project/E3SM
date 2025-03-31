set(ACE_SRCS
  eamxx_ace_process_interface.cpp
  ace_grids_manager.cpp
)

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
  add_library(ace ${ACE_SRCS})
  target_compile_definitions(ace PUBLIC EAMXX_HAS_ACE)
  target_include_directories(ace SYSTEM PUBLIC ${TORCH_INCLUDE_DIRS})

  target_link_libraries(ace PUBLIC
    physics_share
    scream_share
    ${TORCH_LIBRARIES}
  )

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
