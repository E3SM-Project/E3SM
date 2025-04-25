set(INFRA_SRCS
  gw_test_data.cpp
  gw_iso_c.f90
)

add_library(gw_test_infra ${INFRA_SRCS})
target_link_libraries(gw_test_infra gw)
target_include_directories(gw_test_infra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
