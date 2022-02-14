set(DIAGNOSTIC_SRCS
 potential_temperature.cpp
)

set(DIAGNOSTIC_HEADERS
  potential_temperature.hpp
) 

add_library(diagnostics ${DIAGNOSTIC_SRCS})
target_include_directories(diagnostics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(diagnostics scream_share)
target_compile_options(diagnostics PUBLIC)

add_subdirectory(tests)
