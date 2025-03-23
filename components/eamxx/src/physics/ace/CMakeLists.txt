set(ACE_SRCS
  eamxx_ace_process_interface.cpp
)

set(CMAKE_PREFIX_PATH "/global/u1/m/mahf708/.micromamba/envs/fme/lib/python3.10/site-packages/torch/share/cmake/Torch")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
find_package(Torch REQUIRED)

set(ACE_LIBS "ace")
add_library(ace ${ACE_SRCS})
target_compile_definitions(ace PUBLIC EAMXX_HAS_ACE)
target_include_directories(ace PUBLIC ${TORCH_INCLUDE_DIRS})
target_link_libraries(ace PUBLIC physics_share scream_share "${TORCH_LIBRARIES}")
target_link_libraries(ace PUBLIC scream)
target_link_options(ace PUBLIC "-Wl,--copy-dt-needed-entries")
target_link_libraries(eamxx_physics INTERFACE ace)
