set(PHYSICS_SHARE_SRCS
  physics_share_f2c.F90
  physics_share.cpp
  physics_test_data.cpp
  eamxx_trcmix.cpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU)
  list(APPEND PHYSICS_SHARE_SRCS
  physics_saturation.cpp)
endif()

add_library(physics_share ${PHYSICS_SHARE_SRCS})
set_target_properties(physics_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(physics_share PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_link_libraries(physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE physics_share)
endif()
