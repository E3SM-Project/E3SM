add_library(spa eamxx_spa.cpp)
target_link_libraries(spa physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
