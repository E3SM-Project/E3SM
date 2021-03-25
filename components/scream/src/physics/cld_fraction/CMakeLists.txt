set(CLDFRAC_SRCS
  atmosphere_cld_fraction.cpp
  cld_fraction.cpp 
)

set(CLDFRAC_HEADERS
  atmosphere_cld_fraction.hpp
  cld_fraction_functions.hpp
  cld_fraction_main_impl.hpp
)

add_library(cldfrac ${CLDFRAC_SRCS})
target_include_directories(cldfrac PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(cldfrac physics_share scream_share)
target_compile_options(cldfrac PUBLIC)

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()
