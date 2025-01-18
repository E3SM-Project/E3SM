add_library(eamxx_physics INTERFACE)

# Build shared physics functions first, so other parametrizations can use them
add_subdirectory(share)

# Add all physics parametrizations subfolders
add_subdirectory(p3)
if (SCREAM_DOUBLE_PRECISION)
  add_subdirectory(rrtmgp)
  add_subdirectory(cosp)
  add_subdirectory(tms)
  add_subdirectory(iop_forcing)
else()
  message(STATUS "WARNING: RRTMGP, COSP, TMS, and IOPForcing only supported for double precision builds; skipping")
endif()
add_subdirectory(shoc)
add_subdirectory(cld_fraction)
if (SCREAM_ENABLE_ML_CORRECTION)
  add_subdirectory(ml_correction)
endif()
add_subdirectory(spa)
add_subdirectory(nudging)
if (SCREAM_ENABLE_MAM)
  add_subdirectory(mam)
endif()

