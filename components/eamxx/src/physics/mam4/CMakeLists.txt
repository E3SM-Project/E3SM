add_library(mam4 atmosphere_aerosols.cpp)
add_dependencies(mam4 mam4xx_proj)
target_include_directories(mam4 PRIVATE
  ${PROJECT_BINARY_DIR}/externals/haero/include
  ${PROJECT_BINARY_DIR}/externals/mam4xx/include
)
target_link_libraries(mam4 physics_share scream_share mam4xx)

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()
