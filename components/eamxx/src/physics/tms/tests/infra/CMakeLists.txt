set(INFRA_SRCS
  tms_test_data.cpp
)

add_library(tms_test_infra ${INFRA_SRCS})
target_link_libraries(tms_test_infra PUBLIC tms scream_test_support)
target_include_directories(tms_test_infra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
