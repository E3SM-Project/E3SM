# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT ${SCREAM_BASELINES_ONLY})
  include (ScreamUtils)

  # Test control folder
  CreateUnitTest(ping_pong "ping_pong_test.cpp" "scream_control;scream_share" LABELS "driver")

  # Copy yaml input file to run directory
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ping_pong.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/ping_pong.yaml COPYONLY)

  # Unit test the ad
  CreateUnitTest(ad_ut "ad_tests.cpp" "scream_control;scream_share" LABELS "driver")

  # Copy yaml input file to run directory
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ad_tests.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/ad_tests.yaml COPYONLY)
endif()
