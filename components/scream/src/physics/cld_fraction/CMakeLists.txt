set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

set(CLDFRAC_SRCS
  atmosphere_cldfraction.cpp
  cld_fraction.cpp 
)

set(CLDFRAC_HEADERS
  atmosphere_cldfraction.hpp
  cldfraction_functions.hpp
  cldfraction_main_impl.hpp
)

add_library(cldfrac ${CLDFRAC_SRCS})
target_include_directories(cldfrac PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_include_directories(cldfrac SYSTEM PUBLIC ${CIMEROOT}/src/share/include)
target_link_libraries(cldfrac physics_share scream_share)
target_compile_options(cldfrac PUBLIC)

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()
