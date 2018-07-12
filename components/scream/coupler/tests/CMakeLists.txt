
set(COUPLER_TEST_SRCS
  coupler_tests.cpp
)

add_executable(coupler_tests ${COUPLER_TEST_SRCS})

target_include_directories(coupler_tests PUBLIC ${CATCH_INCL_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(coupler_tests coupler)

add_test(run_coupler_tests coupler_tests)
