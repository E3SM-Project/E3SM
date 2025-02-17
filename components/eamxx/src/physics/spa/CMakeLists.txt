add_library(spa eamxx_spa_process_interface.cpp)
target_compile_definitions(spa PUBLIC EAMXX_HAS_SPA)
target_link_libraries(spa physics_share scream_share)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE spa)
endif()
