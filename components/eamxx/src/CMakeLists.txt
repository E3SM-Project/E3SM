add_subdirectory(share)

add_subdirectory(dynamics)
add_subdirectory(physics)
add_subdirectory(diagnostics)
add_subdirectory(control)
if (PROJECT_NAME STREQUAL "E3SM")
  add_subdirectory(mct_coupling)
endif()

option (EAMXX_ENABLE_PYSCREAM "Whether to enable python interface to eamxx" OFF)
if (EAMXX_ENABLE_PYSCREAM)
  add_subdirectory(python)
endif()
