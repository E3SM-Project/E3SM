set(CONTROL_TEST_SOURCES
  control_tests.cpp
)

link_directories(${SCREAM_LIBRARY_DIRS} ${SCREAM_TPL_LIBRARY_DIRS})

add_executable(control_tests ${CONTROL_TEST_SOURCES})
target_include_directories(control_tests PUBLIC ${CATCH_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(control_tests scream_control scream_share ${SCREAM_TPL_LIBRARIES})

add_test(run_control_tests control_tests)
