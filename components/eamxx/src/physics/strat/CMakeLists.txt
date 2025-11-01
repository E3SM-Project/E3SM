add_library(strat_chem eamxx_linoz_process_interface.cpp)
target_compile_definitions(strat_chem PUBLIC EAMXX_HAS_STRAT_CHEM)
target_link_libraries(strat_chem eamxx_physics_share scream_share mam)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE strat_chem)
endif()
