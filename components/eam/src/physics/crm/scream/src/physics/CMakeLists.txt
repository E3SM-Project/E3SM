# Build shared physics functions first, so other parametrizations can use them
add_subdirectory(share)

# Add all physics parametrizations subfolders
add_subdirectory(p3)
if (SCREAM_DOUBLE_PRECISION)
  add_subdirectory(rrtmgp)
else()
  message(STATUS "WARNING: RRTMGP only supported for double precision builds; skipping")
endif()
add_subdirectory(shoc)
if (SCREAM_ENABLE_EXPERIMENTAL)
  add_subdirectory(zm)
endif()
add_subdirectory(cld_fraction)
add_subdirectory(spa)
