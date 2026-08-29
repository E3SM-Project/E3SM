cmake_minimum_required(VERSION 3.3)
cmake_policy(SET CMP0057 NEW)

project(SCRIPTS_TESTS)

enable_testing()
include(CTest)

add_test(NAME scripts_tests
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../scripts-tests -f -m $ENV{SCREAM_MACHINE})

add_test(NAME cime-nml-tests
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../cime-nml-tests)
