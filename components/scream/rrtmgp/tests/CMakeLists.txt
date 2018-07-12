
set(RRTMGP_TEST_SRCS
  rrtmgp_tests.cpp
)

add_executable(rrtmgp_tests ${RRTMGP_TEST_SRCS})

target_include_directories(rrtmgp_tests PUBLIC ${CATCH_INCL_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(rrtmgp_tests rrtmgp)

add_test(run_rrtmgp_tests rrtmgp_tests)
