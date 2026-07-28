cmake_minimum_required(VERSION 3.16)

enable_testing()

set(COUPLER_API_TEST_SOURCES
  test_coupler_api.F90
)

add_executable(test_coupler_api ${COUPLER_API_TEST_SOURCES})

target_link_libraries(test_coupler_api PRIVATE
  E3SM::coupler_api
)
add_test(NAME coupler_api_basic
         COMMAND test_coupler_api)
