# Here, we process a few shared libraries shipped with e3sm.
# For each of them, we either
#   - build the library (standalone build) or
#   - wrap pre-built library in a CMake target (CIME build)

# First pioc/piof, since we link against it in csm_share (at least in CIME build)
include (${SCREAM_BASE_DIR}/cmake/tpls/Scorpio.cmake)
CreateScorpioTargets()

# Then csm_share
include (${SCREAM_BASE_DIR}/cmake/tpls/CsmShare.cmake)
CreateCsmShareTarget()

if (SCREAM_CIME_BUILD)
  # For CIME runs, wrap mct in a target too
  include (${SCREAM_BASE_DIR}/cmake/tpls/Mct.cmake)
  CreateMctTarget()
endif()
