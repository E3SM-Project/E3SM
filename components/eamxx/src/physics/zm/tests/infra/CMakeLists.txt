set(INFRA_SRCS
  # zm_test_data.cpp
  zm_iso_c.f90
)

add_library(zm_test_infra ${INFRA_SRCS})
target_link_libraries(zm_test_infra zm)
target_include_directories(zm_test_infra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
