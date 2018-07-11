
file(GLOB SCREAM_TEST_SRCS "*.cpp")

add_executable(scream_tests ${SCREAM_TEST_SRCS})

target_include_directories(scream_tests PUBLIC ${CATCH_INCL_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(scream_tests coupler p3 rrtmgp shoc)

add_test(run_scream_tests scream_tests)
