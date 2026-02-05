add_library(spc eamxx_spc_process_interface.cpp)
target_compile_definitions(spc PUBLIC EAMXX_HAS_SPC)
target_link_libraries(spc eamxx_physics_share scream_share)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE spc)
endif()
