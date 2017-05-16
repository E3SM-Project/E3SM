cmake_minimum_required(VERSION 2.8)

include(CTest)

execute_process(COMMAND "python" "list_tests" OUTPUT_VARIABLE STR_TESTS
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_STRIP_TRAILING_WHITESPACE)

separate_arguments(TEST_LIST UNIX_COMMAND ${STR_TESTS})

foreach(ATEST ${TEST_LIST})
    add_test(NAME ${ATEST} COMMAND "python" "-m" "unittest" "-v" ${ATEST}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach(ATEST)
