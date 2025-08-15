add_library(iop_forcing eamxx_iop_forcing_process_interface.cpp)
target_compile_definitions(iop_forcing PUBLIC EAMXX_HAS_IOP_FORCING)
target_link_libraries(iop_forcing physics_share scream_share)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE iop_forcing)
endif()
