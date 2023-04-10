cmake_minimum_required(VERSION 3.3)
cmake_policy(SET CMP0057 NEW)

project(SCRIPTS_TESTS)

enable_testing()
include(CTest)

add_test(NAME scripts_tests_tas COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../gather-all-data -l -m $ENV{SCREAM_MACHINE} "${CMAKE_CURRENT_SOURCE_DIR}/../scripts-tests -f -m \$machine TestTestAllScream")
add_test(NAME scripts_tests_jenkins COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../scripts-tests -j -m $ENV{SCREAM_MACHINE} TestTestAllScream)
