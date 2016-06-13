cmake_minimum_required(VERSION 2.8)

include(CTest)

add_test (NAME cime_regression_default COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts_regression_tests)
add_test (NAME cime_regression_no_batch COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts_regression_tests --no-batch)
