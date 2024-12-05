set(INFRA_SRCS
  p3_data.cpp
  p3_ic_cases.cpp
  p3_main_wrap.cpp
  p3_test_data.cpp
)

#crusher change
if (Kokkos_ENABLE_HIP)
set_source_files_properties(p3_test_data.cpp  PROPERTIES COMPILE_FLAGS -O0)
endif()

add_library(p3_test_infra ${INFRA_SRCS})
target_link_libraries(p3_test_infra p3)
target_include_directories(p3_test_infra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
