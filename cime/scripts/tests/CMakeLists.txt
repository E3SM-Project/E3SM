#
# Environment variables CIME_COMPILER and CIME_MPILIB
# can be used to send --compiler and --mpilib settings to scripts_regression_tests.py
#
PROJECT(cime)
cmake_minimum_required(VERSION 2.8)

include(CTest)
if (DEFINED ENV{PYTHON})
  set(PYTHON $ENV{PYTHON})
else()
  set(PYTHON "python")
endif()

execute_process(COMMAND ${PYTHON} "--version" OUTPUT_VARIABLE PY_VER
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_STRIP_TRAILING_WHITESPACE)

MESSAGE("Python version is " ${PY_VER})		     
execute_process(COMMAND ${PYTHON} "list_tests" OUTPUT_VARIABLE STR_TESTS
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_STRIP_TRAILING_WHITESPACE)
separate_arguments(TEST_LIST UNIX_COMMAND ${STR_TESTS})

set(args "")
if (DEFINED ENV{CIME_COMPILER})
  set(args "--compiler $ENV{CIME_COMPILER}")
endif()
if (DEFINED ENV{CIME_MPILIB})
  set(args "${args} --mpilib $ENV{CIME_MPILIB}")
endif()
separate_arguments(ARG_LIST UNIX_COMMAND ${args})

foreach(ATEST ${TEST_LIST})
    # this assignment prevents quotes being added to testname in add_test
    set(fulltest "${ATEST}")
    add_test(NAME ${ATEST} COMMAND ./scripts_regression_tests.py -v ${fulltest} ${ARG_LIST}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach(ATEST)
