add_library(spa eamxx_spa_process_interface.cpp)
target_compile_definitions(spa PUBLIC EAMXX_HAS_SPA)
target_link_libraries(spa physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
