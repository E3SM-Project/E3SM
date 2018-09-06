set(SCREAM_TEST_SRCS
  ../util/scream_catch_main.cpp
  scream_tests.cpp
)

link_directories(${SCREAM_LIBRARY_DIRS})
add_executable(scream_tests ${SCREAM_TEST_SRCS})

target_include_directories(scream_tests PUBLIC ${CATCH_INCL_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(scream_tests coupler p3 rrtmgp shoc)
set_target_properties(scream_tests PROPERTIES LINK_FLAGS "${SCREAM_LINK_FLAGS}")
target_link_libraries(scream_tests scream_util ${SCREAM_LIBRARIES})

add_test(run_scream_tests scream_tests)
