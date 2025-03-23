set(ACE_SRCS
  eamxx_ace_process_interface.cpp
)

# Use a configurable path for PyTorch with a default value
if(NOT DEFINED TORCH_CMAKE_PATH)
  set(TORCH_CMAKE_PATH "/global/u1/m/mahf708/.micromamba/envs/fme/lib/python3.10/site-packages/torch/share/cmake/Torch" 
      CACHE PATH "Path to PyTorch CMake files")
endif()
set(CMAKE_PREFIX_PATH ${TORCH_CMAKE_PATH})
find_package(Torch REQUIRED)

include(ScreamUtils)
add_library(ace ${ACE_SRCS})
target_compile_definitions(ace PUBLIC EAMXX_HAS_ACE)
target_include_directories(ace PUBLIC ${TORCH_INCLUDE_DIRS})

target_link_libraries(ace PUBLIC
  physics_share
  scream_share
  ${TORCH_LIBRARIES}
)

if (TARGET eamxx_physics)
  target_link_libraries(eamxx_physics INTERFACE ace)
endif()
