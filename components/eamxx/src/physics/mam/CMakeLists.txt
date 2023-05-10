add_library(mam eamxx_mam_microphysics.cpp)
add_dependencies(mam mam4xx_proj)
target_include_directories(mam PRIVATE
  ${PROJECT_BINARY_DIR}/externals/haero/include
  ${PROJECT_BINARY_DIR}/externals/mam4xx/include
)
target_link_libraries(mam physics_share scream_share mam4xx)

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()
