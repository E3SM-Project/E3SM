set(TMS_SRCS
  eamxx_tms_process_interface.cpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND TMS_SRCS
    eti/compute_tms.cpp
  )
endif()

# Create tms library
add_library(tms ${TMS_SRCS})
target_link_libraries(tms physics_share scream_share)
target_compile_definitions(tms PUBLIC EAMXX_HAS_TMS)
target_include_directories(tms PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/impl
)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE tms)
endif()
