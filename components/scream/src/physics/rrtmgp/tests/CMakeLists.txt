include(ScreamUtils)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  CreateUnitTest(rrtmgp_tests rrtmgp_tests.cpp "rrtmgp;scream_share")
endif()
