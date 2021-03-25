set(CLDFRAC_SRCS
  atmosphere_cld_fraction.cpp
  cld_fraction.cpp 
)

set(CLDFRAC_HEADERS
  atmosphere_cld_fraction.hpp
  cld_fraction_functions.hpp
  cld_fraction_main_impl.hpp
)

add_library(cld_fraction ${CLDFRAC_SRCS})
target_include_directories(cld_fraction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_link_libraries(cld_fraction physics_share scream_share)
target_compile_options(cld_fraction PUBLIC)

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()
