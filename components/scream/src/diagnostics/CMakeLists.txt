set(DIAGNOSTIC_SRCS
  potential_temperature.cpp
  atm_density.cpp
  exner.cpp
  virtual_temperature.cpp
)

add_library(diagnostics ${DIAGNOSTIC_SRCS})
target_include_directories(diagnostics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(diagnostics PUBLIC scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
