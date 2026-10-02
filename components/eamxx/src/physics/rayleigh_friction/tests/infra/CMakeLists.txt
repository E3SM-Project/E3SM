set(INFRA_SRCS
  rayleigh_friction_test_data.cpp
)

add_library(rayleigh_friction_test_infra ${INFRA_SRCS})
target_link_libraries(rayleigh_friction_test_infra PUBLIC rayleigh_friction eamxx_test_support)
target_include_directories(rayleigh_friction_test_infra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
