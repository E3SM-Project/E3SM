set(SCREAM_CONTROL_SOURCES
  atmosphere_driver.cpp
  surface_coupling.cpp
)

set(SCREAM_CONTROL_HEADERS
  atmosphere_driver.hpp
  surface_coupling.hpp
)

link_directories(${SCREAM_TPL_LIBRARY_DIRS})

add_library(scream_control ${SCREAM_CONTROL_SOURCES})
target_include_directories(scream_control PUBLIC ${SCREAM_SRC_DIR})
target_link_libraries(scream_control scream_share ${SCREAM_TPL_LIBRARIES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
