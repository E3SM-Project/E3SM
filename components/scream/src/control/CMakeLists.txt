set(SCREAM_CONTROL_SOURCES
  atmosphere_driver.cpp
)

set(SCREAM_CONTROL_HEADERS
  atmosphere_driver.hpp
)

add_library(scream_control ${SCREAM_CONTROL_SOURCES})
target_link_libraries(scream_control scream_share scream_io)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
