add_subdirectory(share)

# Also handle csm_share. If this is a CIME build,
# it will wrap the lib in a cmake target, otherwise
# it will build a small lib (same name in both cases).
include(${SCREAM_TPLS_MODULE_DIR}/CsmShare.cmake)
CreateCsmShareTarget()

add_subdirectory(dynamics)
add_subdirectory(physics)
add_subdirectory(diagnostics)
add_subdirectory(control)
if (PROJECT_NAME STREQUAL "E3SM")
  add_subdirectory(mct_coupling)
endif()
