if (NOT SCREAM_LIB_ONLY)
  # Create test_support lib. Subfolders will add sources as we scan them
  add_library(eamxx_test_support)
  target_link_libraries (eamxx_test_support PUBLIC Catch2::Catch2)
endif()

add_subdirectory(share)

add_subdirectory(dynamics)
add_subdirectory(physics)
add_subdirectory(control)
if (PROJECT_NAME STREQUAL "E3SM")
  add_subdirectory(mct_coupling)
endif()

option (EAMXX_ENABLE_PYSCREAM "Whether to enable interfaces to call eamxx from python" OFF)
if (EAMXX_ENABLE_PYSCREAM)
  add_subdirectory(python)
endif()
