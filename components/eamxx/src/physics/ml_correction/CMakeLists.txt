set(MLCORRECTION_SRCS
  eamxx_ml_correction_process_interface.cpp
)

set(MLCORRECTION_HEADERS
  eamxx_ml_correction_process_interface.hpp
)

add_library(ml_correction ${MLCORRECTION_SRCS})
target_include_directories(ml_correction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(ml_correction physics_share scream_share)
target_compile_options(ml_correction PUBLIC)
