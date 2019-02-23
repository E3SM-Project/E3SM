set(SCREAM_CONTROL_SOURCES
  cxx_f90_interface.cpp
  atmosphere_driver.cpp
  atmosphere_process_factory.cpp
  atmosphere_process_group.cpp
  surface_coupling.cpp
)

set(SCREAM_CONTROL_HEADERS
  atmosphere_driver.hpp
  atmosphere_process_factory.hpp
  atmosphere_process_group.hpp
  surface_coupling.hpp
)

link_directories(${SCREAM_LIBRARY_DIRS} ${SCREAM_TPL_LIBRARY_DIRS})

add_library(scream_control ${SCREAM_CONTROL_SOURCES})
target_include_directories(scream_control PUBLIC ${SCREAM_SRC_DIR})
target_link_libraries(scream_control scream_dynamics scream_share ${SCREAM_TPL_LIBRARIES})

add_subdirectory(tests)
