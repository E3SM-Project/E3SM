set(SCREAM_TEST_SRCS
  ${SCREAM_SRC_DIR}/share/util/scream_catch_main.cpp
  scream_tests.cpp
)

link_directories(${SCREAM_LIBRARY_DIRS} ${SCREAM_TPL_LIBRARY_DIRS})

add_executable(scream_tests ${SCREAM_TEST_SRCS})
target_include_directories(scream_tests PUBLIC ${SCREAM_INCLUDE_DIRS} ${SCREAM_TPL_INCLUDE_DIRS} ${CATCH_INCLUDE_DIR})
target_link_libraries(scream_tests scream_control p3 rrtmgp shoc scream_share ${SCREAM_TPL_LIBRARIES})
set_target_properties(scream_tests PROPERTIES LINK_FLAGS "${SCREAM_LINK_FLAGS}")

add_test(run_scream_tests scream_tests)
