add_library (eamxx_physics_share
  eamxx_trcmix.cpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU)
  target_sources(eamxx_physics_share PRIVATE
    physics_saturation.cpp)
endif()

set_target_properties(eamxx_physics_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(eamxx_physics_share PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_link_libraries(eamxx_physics_share PUBLIC
  eamxx_core
  eamxx_utils
  ekat::Pack)

if (NOT SCREAM_LIB_ONLY)
  # Add source files that are only used in testing
  target_sources(eamxx_physics_share PRIVATE
    physics_share_f2c.F90
    physics_share.cpp
  )

  # This file does not add any other dependency to eamxx_test_support,
  # so it makes sense to add it to that lib (unlike mesh_free_grids_manager.cpp,
  # which require linking against lots of other eamxx libs)
  target_sources(eamxx_test_support PRIVATE
    physics_test_data.cpp
  )

  add_subdirectory(tests)
endif()
