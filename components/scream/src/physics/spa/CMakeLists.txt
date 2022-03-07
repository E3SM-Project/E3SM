add_library(spa atmosphere_prescribed_aerosol.cpp)
target_link_libraries(spa physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
