set(MLNUDGING_SRCS
  atmosphere_ml_nudging.cpp
)

set(MLNUDGING_HEADERS
  atmosphere_ml_nudging.hpp
)

add_library(ml_nudging ${MLNUDGING_SRCS})
target_include_directories(ml_nudging PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(ml_nudging physics_share scream_share)
target_compile_options(ml_nudging PUBLIC)
