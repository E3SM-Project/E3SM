add_library(nudging atmosphere_nudging.cpp)
target_include_directories(nudging PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../share)
target_link_libraries(nudging physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
