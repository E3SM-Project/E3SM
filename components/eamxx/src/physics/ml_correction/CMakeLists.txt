set(MLCORRECTION_SRCS
  atmosphere_ml_correction.cpp
)

set(MLCORRECTION_HEADERS
  atmosphere_ml_correction.hpp
)

add_library(ml_correction ${MLCORRECTION_SRCS})
target_include_directories(ml_correction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(ml_correction physics_share scream_share)
target_compile_options(ml_correction PUBLIC)
