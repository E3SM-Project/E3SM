cmake_minimum_required(VERSION 2.8)

include(CTest)

add_test (NAME acme_scripts_regression COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts_regression_tests)
