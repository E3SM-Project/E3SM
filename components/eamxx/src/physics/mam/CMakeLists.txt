add_library(mam
  eamxx_mam_microphysics_process_interface.cpp
  eamxx_mam_optics_process_interface.cpp)
target_compile_definitions(mam PUBLIC EAMXX_HAS_MAM)
add_dependencies(mam mam4xx_proj)
target_include_directories(mam PUBLIC
  ${PROJECT_BINARY_DIR}/externals/haero/include
  ${PROJECT_BINARY_DIR}/externals/mam4xx/include
)
target_link_libraries(mam PUBLIC physics_share scream_share mam4xx haero)

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()
