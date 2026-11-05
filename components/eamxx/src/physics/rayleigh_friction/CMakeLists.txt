set(RAYLEIGH_FRICTION_SRCS
  eamxx_rayleigh_friction_process_interface.cpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND RAYLEIGH_FRICTION_SRCS
    eti/rayleigh_friction.cpp
  )
endif()

# Create rayleigh_friction library
add_library(rayleigh_friction ${RAYLEIGH_FRICTION_SRCS})
target_link_libraries(rayleigh_friction eamxx_physics_share scream_share)
target_compile_definitions(rayleigh_friction PUBLIC EAMXX_HAS_RAYLEIGH_FRICTION)
target_include_directories(rayleigh_friction PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/impl
)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE rayleigh_friction)
endif()
