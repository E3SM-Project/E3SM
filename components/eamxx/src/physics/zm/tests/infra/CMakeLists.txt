set(INFRA_SRCS
  zm_iso_c.f90
  zm_test_data.cpp
)

add_library(zm_test_infra ${INFRA_SRCS})
target_link_libraries(zm_test_infra zm scream_test_support)
target_include_directories(zm_test_infra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
