
file(GLOB P3_TEST_SRCS "*.cpp")

add_executable(p3_tests ${P3_TEST_SRCS})

target_include_directories(p3_tests PUBLIC ${CATCH_INCL_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(p3_tests p3)

add_test(run_p3_tests p3_tests)
