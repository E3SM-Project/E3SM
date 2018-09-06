link_directories(${SCREAM_LIBRARY_DIRS})

add_executable(p3_tests
  ${SCREAM_SRC_DIR}/share/util/scream_catch_main.cpp
  p3_tests.cpp)
add_executable(p3_run_and_cmp
  p3_run_and_cmp.cpp)

add_custom_target(p3_baseline
  COMMAND $<TARGET_FILE:p3_run_and_cmp> -g ${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline)

add_dependencies(baseline p3_baseline)

foreach (exe p3_tests p3_run_and_cmp)
  target_include_directories(${exe} PUBLIC ${CATCH_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
  target_link_libraries(${exe} p3)
  set_target_properties(${exe} PROPERTIES LINK_FLAGS "${SCREAM_LINK_FLAGS}")
  target_link_libraries(${exe} scream_share ${SCREAM_LIBRARIES})
endforeach ()

add_test(run_p3_tests p3_tests)
configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v2.8.2 p3_lookup_table_1.dat-v2.8.2 COPYONLY)

add_test(p3_regression p3_run_and_cmp ${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline)
