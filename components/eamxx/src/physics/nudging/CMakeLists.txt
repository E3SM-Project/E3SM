add_library(nudging eamxx_nudging_process_interface.cpp)
target_compile_definitions(nudging PUBLIC EAMXX_HAS_NUDGING)
target_link_libraries(nudging physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE nudging)
endif()
