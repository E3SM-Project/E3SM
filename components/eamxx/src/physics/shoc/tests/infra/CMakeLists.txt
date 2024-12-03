set(INFRA_SRCS
  shoc_data.cpp
  shoc_ic_cases.cpp
  shoc_main_wrap.cpp
  shoc_test_data.cpp
)

add_library(shoc_test_infra ${INFRA_SRCS})
target_link_libraries(shoc_test_infra shoc)
target_include_directories(shoc_test_infra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
