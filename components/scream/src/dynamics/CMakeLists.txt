set(SCREAM_DYNAMICS_SOURCES
  atmosphere_dynamics.cpp
)

set(SCREAM_DYNAMICS_HEADERS
  atmosphere_dynamics.hpp
)

link_directories(${SCREAM_LIBRARY_DIRS} ${SCREAM_TPL_LIBRARY_DIRS})

add_library(scream_dynamics ${SCREAM_DYNAMICS_SOURCES})
target_include_directories(scream_dynamics PUBLIC ${SCREAM_SRC_DIR})
target_link_libraries(scream_dynamics scream_share ${SCREAM_TPL_LIBRARIES})
