set(INFRA_SRCS
  gw_test_data.cpp
)

add_library(gw_test_infra ${INFRA_SRCS})
target_link_libraries(gw_test_infra PUBLIC gw eamxx_test_support)
target_include_directories(gw_test_infra PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
